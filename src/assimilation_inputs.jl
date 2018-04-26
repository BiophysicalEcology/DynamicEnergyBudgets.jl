using PlantPhysiology

# function photomoduleC3(s::S) where S
#     net, gross, dark = Photosynthesis.photosynthesisC3(s.assim_state, s.params) 
#     return net
# end

function maestra(o)
    o.assim_state.aleaf * o.params.SLA * o.params.w_V
end

function photosynthesis(f::KooijmanPhotosynthesis, o)
    p = o.params
    j1_l = half_saturation(p.J_L_F, p.J_L_K, p.j_L_Amax)
    j1_c = half_saturation(p.X_C, p.K_C, p.j_C_Amax)
    j1_o = half_saturation(p.X_O, p.K_O, p.j_O_Amax)

    # photorespiration.
    bound_o = j1_o/p.k_O_binding
    bound_c = j1_c/p.k_C_binding

    # c flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l) # assimilated c-reserve (mol/time)/(1 + mol/mol + mol/mol + dimless) = mol/time
end

function photosynthesis(f::KooijmanSLAPhotosynthesis, o)
    j1_l = half_saturation(f.J_L_F, f.J_L_K, f.j_L_Amax) # mol.s⁻¹m⁻²
    j1_c = half_saturation(f.X_C, f.K_C, f.j_C_Amax) # mol.s⁻¹m⁻²
    j1_o = half_saturation(f.X_O, f.K_O, f.j_O_Amax) # mol.s⁻¹m⁻²

    mass_area_coef = o.paramo.w_V * f.SLA # m⁻².g * g.mol⁻¹ = m².mol⁻¹
    j1_l = j1_l * mass_area_coef
    j1_c = j1_c * mass_area_coef
    j1_o = j1_o * mass_area_coef # mol.s⁻¹m⁻² * m².mol⁻¹ = mol.mol.s⁻¹

    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol.s⁻¹/molo.s⁻¹ = mol/mol
    bound_c = j1_c/f.k_C_binding # mol.s⁻¹/molo.s⁻¹ = mol/mol

    # c flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co) # mol.m⁻²s⁻¹/mol.m⁻²s⁻¹ - mol.m⁻²s⁻¹/mol.m⁻²s⁻¹ = dimless

    j_c_intake / (1 + bound_c + bound_o + co_l) # assimilated c-reserve (mol/time)/(1 + mol/mol + mol/mol + mol/mol) = mol/time
end

function feedback!(o, f::Autophagy, u::AbstractState, r)
    ap = u[V] * (1 - half_saturation(r, f.K_autophagy, oneunit(r)))
    o.J[:C,:gro] += ap/o.params.y_E_CH_NO
    o.J[:N,:gro] += ap/o.params.y_E_EN
    o.J[:V,:gro] -= ap
    nothing
end
function feedback!(o, f::Autophagy, u::AbstractStateE, r)
    ap = u[V] * (1 - half_saturation(r, f.K_autophagy, oneunit(r)))
    o.J[:E,:gro] += ap
    o.J[:V,:gro] -= ap
    nothing
end

"""
    area_mass_kooijman(uV, Vref, Vscaling)
Area mass scaling for plants.
(uV / Vref)^(1/3 - (-uV / Vscaling)^b) DEB book p. 134
"""
function area_mass_kooijman(uV, Vref, Vscaling)
    uV > zero(uV) || return zero(uV)
    (uV / Vref)^(-uV / Vscaling)
    # (uV / Vref)^(1/3 - (-uV / Vscaling)^b) DEB book p. 134
end

function area_mass_linear(uV, Vref, Vscaling)
    return uV
end


function assimilation!(o, on, f::CarbonAssimilation, u::AbstractStateCN, p)::Void
    p = o.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:N,:rej] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1])
        return nothing
    end
    J1_EC_ass = photosynthesize(f.formulation, o) * u.V * o.A

    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * sn.J1[:N,:rej] 

    (o.J[:C,:ass], o.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)
    return nothing
end
function assimilation!(o, on, f::CarbonAssimilation, u::AbstractStateCNE)::Void
    p = o.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:N,:rej] = o.J[:E,:ass] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1])
        return nothing
    end

    J1_EC_ass = photosynthesize(f.formulation, o) * u.V * o.A

    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * sn.J1[:N,:rej] 

    # Merge rejected N from root and photosynthesized C into reserves
    (o.J[:C,:ass], o.J[:N,:ass], o.J[:E,:ass]) = 
        synthesizing_unit(J1_EC_ass, J1_EN_ass, p.y_E_CH_NO, p.y_E_EN)
    return nothing
end
function assimilation!(o, on, f::NH4_NO3_Assimilation, u::AbstractStateCN)::Void
    p = o.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:C,:rej] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1]) 
        return nothing
    end
    # Arriving nitrogen
    (J_N_ass, J_NO_ass, J_NH_ass) = uptake_nitrogen(f.formulation, o, on)

    # Rejected C-reserve from shoot
    J1_EC_ass = -p.y_EC_ECT * sn.J1[:C,:rej]

    θNH = J_NH_ass/J_N_ass                           # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                    # Fraction of nitrate in arriving N-flux
    y_E_CH = θNH * p.y_E_CH_NH + θNO * p.y_E_CH_NO   # Yield coefficient from C-reserve to reserve

    (o.J[:C,:ass], o.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    # How to deal with this without general reserve?
    # o.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * o.J[:E,:ass]) * 1/p.n_N_EN
    return nothing
end
function assimilation!(o, on, f::NH4_NO3_Assimilation, u::AbstractStateCNE)
    p = o.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:C,:rej] = o.J[:E,:ass] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1]) 
        return nothing
    end
    # Arriving nitrogen
    (J_N_ass, J_NO_ass, J_NH_ass) = uptake_nitrogen(f.formulation, o, on)

    J1_EC_ass = -p.y_EC_ECT * sn.J1[:C,:rej] # Rejected C-reserve from shoot

    θNH = J_NH_ass/J_N_ass                          # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                   # Fraction of nitrate in arriving N-flux
    y_E_CH = θNH * p.y_E_CH_NH + θNO * p.y_E_CH_NO  # Yield coefficient from C-reserve to reserve

    # Merge rejected C from shoot and uptaken N into reserves
    (o.J[:C,:ass], o.J[:N,:ass], o.J[:E,:ass]) = 
        synthesizing_unit(J1_EC_ass, J_N_ass, y_E_CH, 1/p.n_N_E)

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    o.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * o.J[:E,:ass]) * 1/p.n_N_EN
    return nothing
end

function uptake_nitrogen(f::Kooijman_NH4_NO3_Assimilation, o, on)
    p = o.params
    AX_SH = p.X_H * sn.A # mol.L⁻¹ Why A_S? Evapotranspiration based on shoot area, not root area. See p. 199 of DEB book.
    AK_RH = p.K_H * o.A # mol.L⁻¹ TODO: only partial root area is not involved in N uptake [@robinson1991what].
    K1_NH = half_saturation(AK_RH, AX_SH, p.K_NH) # Ammonia saturation mol.mol⁻¹
    K1_NO = half_saturation(AK_RH, AX_SH, p.K_NO) # Nitrate saturation mol.mol⁻¹
    j1_NH = half_saturation(p.X_NH, K1_NH, p.j_NH_Amax) # Arriving ammonia mols.mol⁻¹.s⁻¹
    j1_NO = half_saturation(p.X_NO, K1_NO, p.j_NO_Amax) # Arriving nitrate mols.mol⁻¹.s⁻¹

    # FIXME: why does it need the * mol^-1 ????
    area = o.u.V * o.A * u"mol^-1"  
    J1_NH_ass = area * j1_NH # mol.s⁻¹ = mols.mol⁻¹.s⁻¹
    J_NO_ass = area * j1_NO
    J_N_ass = J1_NH_ass + p.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

function uptake_nitrogen(f::KooijmanSLA_NH4_NO3_Assimilation, o, on)
    p = o.params
    SAX_H = p.X_H * sn.A # mol.L⁻¹ Why A_S? Evapotranspiration based on shoot area, not root area. See p. 199 of DEB book.
    RAK_H = p.K_H * o.A # mol.L⁻¹ TODO: only partial root area is not involved in N uptake [@robinson1991what].
    K1_NH = half_saturation(RAK_H, SAX_H, p.K_NH) # mol.L⁻¹ Ammonia saturation
    K1_NO = half_saturation(RAK_H, SAX_H, p.K_NO) # mol.L⁻¹ Nitrate saturation
    j1_NH = half_saturation(p.X_NH, K1_NH, p.j_NH_Amax) # Arriving ammonia mols.m⁻².s⁻¹ 
    j1_NO = half_saturation(p.X_NO, K1_NO, p.j_NO_Amax) # Arriving nitrate mols.m⁻².s⁻¹

    mass_area_coef = p.w_V * p.SLA # m².mol⁻¹ = m⁻².g * g.mol⁻¹
    area = o.u.V * o.A * mass_area_coef # m² = mol * dimless * m².mol⁻¹
    J1_NH_ass = area * j1_NH # mol.s⁻¹ = m² * mols.m⁻².s⁻¹
    J_NO_ass = area * j1_NO # mol.s⁻¹ = m² * mols.m⁻².s⁻¹
    J_N_ass = J1_NH_ass + p.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

"""
Check if germination has happened. Independent for each organ,
although this may not make sense.
"""
function germinated(M_V, M_Vgerm)::Bool 
    M_V > M_Vgerm 
end


