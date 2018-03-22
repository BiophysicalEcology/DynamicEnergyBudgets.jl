using Photosynthesis

export area_mass_kooijman, shoot_assimilation!, 
photosynthesis_sla, photosynthesis_kooijman, photomoduleC3,
root_assimilation!, nitrogen_uptake, nitrogen_uptake_sla

function photomoduleC3(s::S) where S
    net, gross, dark = Photosynthesis.photosynthesisC3(s.assim_state, s.params) 
    return net
end

function photosynthesis_kooijman(s::S)::Float64 where S
    p = s.params
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

    j1_EC_ass = j_c_intake / (1 + bound_c + bound_o + co_l) # assimilated c-reserve (mol/time)/(1 + mol/mol + mol/mol + dimless) = mol/time
    return j1_EC_ass
end

function photosynthesis_sla(s::S)::Float64 where S
    p = s.params
    j1_l = half_saturation(p.J_L_F, p.J_L_K, p.j_L_Amax) # mol.s⁻¹m⁻²
    j1_c = half_saturation(p.X_C, p.K_C, p.j_C_Amax) # mol.s⁻¹m⁻²
    j1_o = half_saturation(p.X_O, p.K_O, p.j_O_Amax) # mol.s⁻¹m⁻²

    mass_area_coef = p.w_V * p.SLA # m⁻².g * g.mol⁻¹ = m².mol⁻¹
    j1_l = j1_l * mass_area_coef
    j1_c = j1_c * mass_area_coef
    j1_o = j1_o * mass_area_coef # mol.s⁻¹m⁻² * m².mol⁻¹ = mol.mol.s⁻¹

    # photorespiration.
    bound_o = j1_o/p.k_O_binding # mol.s⁻¹/mols.s⁻¹ = mol/mol
    bound_c = j1_c/p.k_C_binding # mol.s⁻¹/mols.s⁻¹ = mol/mol

    # c flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co) # mol.m⁻²s⁻¹/mol.m⁻²s⁻¹ - mol.m⁻²s⁻¹/mol.m⁻²s⁻¹ = dimless

    j1_ec_ass = j_c_intake / (1 + bound_c + bound_o + co_l) # assimilated c-reserve (mol/time)/(1 + mol/mol + mol/mol + mol/mol) = mol/time
    return j1_ec_ass
end

function autophagy!(s::S, u::AbstractState, r::Float64)::Void where S
    ap = u[V] * (1 - half_saturation(r, s.params.K_autophagy, 1.0))
    s.J[:C,:gro] += ap/s.params.y_E_CH_NO * 0.9
    s.J[:N,:gro] += ap/s.params.y_E_EN * 0.9
    s.J[:V,:gro] -= ap
    return nothing
end

function area_mass_kooijman(uV::Float64, Vref::Float64, Vscaling::Float64)::Float64
    uV > 0 || return 0
    (uV / Vref)^(-uV / Vscaling)
    # (uV / Vref)^(1/3 - (-uV / Vscaling)^b) DEB book p. 134
end

function area_mass_linear(uV::Float64, Vref::Float64, Vscaling::Float64)::Float64
    return uV
end

#############################################################################
function shoot_assimilation!(s::S, sn::SN, u::AbstractStateCN, p::T)::Void where {S,SN,T}
    p = s.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:N,:rej] = s.J[:C,:ass] = s.J[:N,:ass] = 0.0
        return nothing
    end
    carbon_per_second = s.functions.assim_sub(s) * u.V * s.A
    J1_EC_ass = carbon_per_second * (60 * 60 * 24 * settings.timestep_days)

    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * sn.J[:N,:rej] 

    (s.J[:C,:ass], s.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)
    return nothing
end

function shoot_assimilation!(s::S, sn::SN, u::AbstractStateCNE, settings::T)::Void where {S,SN,T}
    p = s.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:N,:rej] = s.J[:E,:ass] = s.J[:C,:ass] = s.J[:N,:ass] = 0.0
        return nothing
    end
    carbon_per_second = s.functions.assim_sub(s)
    J1_EC_ass = carbon_per_second * (60 * 60 * 24 * settings.timestep_days)

    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * sn.J[:N,:rej] 

    # Merge rejected N from root and photosynthesized C into reserves
    (s.J[:C,:ass], s.J[:N,:ass], s.J[:E,:ass]) = 
    synthesizing_unit(J1_EC_ass, J1_EN_ass, p.y_E_CH_NO, p.y_E_EN)
    return nothing
end

###############################################################################
function root_assimilation!(s::S, sn::SN, u::AbstractStateCN, settings::T)::Void where {S,SN,T}
    p = s.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:C,:rej] = s.J[:C,:ass] = s.J[:N,:ass] = 0.0 
        return nothing
    end
    # Arriving nitrogen
    uptake_per_second = s.functions.assim_sub(s)
    (J_N_ass, J_NO_ass, J_NH_ass) = uptake_per_second .* (60 * 60 * 24 * settings.timestep_days)
    # Rejected C-reserve from shoot
    J1_EC_ass = -p.y_EC_ECT * sn.J[:C,:rej]

    θNH = J_NH_ass/J_N_ass                                       # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                                # Fraction of nitrate in arriving N-flux
    y_E_CH = θNH * p.y_E_CH_NH + θNO * p.y_E_CH_NO   # Yield coefficient from C-reserve to reserve

    (s.J[:C,:ass], s.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    # How to deal with this without general reserve?
    # s.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * s.J[:E,:ass]) * 1/p.n_N_EN
    return nothing
end

function root_assimilation!(s::S, sn::SN, u::AbstractStateCNE, settings::T)::Void where {S,SN,T}
    p = s.params
    if !germinated(u.V, p.M_Vgerm) 
        sn.J[:C,:rej] = s.J[:E,:ass] = s.J[:C,:ass] = s.J[:N,:ass] = 0.0 
        return nothing
    end
    # Arriving nitrogen
    uptake_per_second = s.functions.assim_sub(s, sn)
    (J_N_ass, J_NO_ass, J_NH_ass) = uptake_per_second .* (60 * 60 * 24 * settings.timestep_days)
    # Rejected C-reserve from shoot
    J1_EC_ass = -p.y_EC_ECT * sn.J[:C,:rej]

    θNH = J_NH_ass/J_N_ass                                       # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                                # Fraction of nitrate in arriving N-flux
    y_E_CH = θNH * p.y_E_CH_NH + θNO * p.y_E_CH_NO   # Yield coefficient from C-reserve to reserve

    # Merge rejected C from shoot and uptaken N into reserves
    (s.J[:C,:ass], s.J[:N,:ass], s.J[:E,:ass]) = 
    synthesizing_unit(J1_EC_ass, J_N_ass, y_E_CH, 1/p.n_N_E)

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    s.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * s.J[:E,:ass]) * 1/p.n_N_EN
    return nothing
end

function nitrogen_uptake(s::S, sn::SN)::NTuple{3,Float64} where {S,SN}
    p = s.params
    AX_SH = p.X_H * sn.A # mol.L⁻¹ Why A_S? Evapotranspiration based on shoot area, not root area. See p. 199 of DEB book.
    AK_RH = p.K_H * s.A # mol.L⁻¹ TODO: only partial root area is not involved in N uptake [@robinson1991what].
    K1_NH = half_saturation(AK_RH, AX_SH, p.K_NH) # Ammonia saturation mol.mol⁻¹
    K1_NO = half_saturation(AK_RH, AX_SH, p.K_NO) # Nitrate saturation mol.mol⁻¹
    j1_NH = half_saturation(p.X_NH, K1_NH, p.j_NH_Amax) # Arriving ammonia mols.mol⁻¹.s⁻¹
    j1_NO = half_saturation(p.X_NO, K1_NO, p.j_NO_Amax) # Arriving nitrate mols.mol⁻¹.s⁻¹

    # TODO: RMA or LMA as a proxy to differentiate root surface area by species
    area = s.u.V * s.A 
    J1_NH_ass = area * j1_NH # mol.s⁻¹ = mols.mol⁻¹.s⁻¹
    J_NO_ass = area * j1_NO
    J_N_ass = J1_NH_ass + p.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

function nitrogen_uptake_sla(s::S, sn::SN)::NTuple{3,Float64} where {S,SN}
    p = s.params
    SAX_H = p.X_H * sn.A # mol.L⁻¹ Why A_S? Evapotranspiration based on shoot area, not root area. See p. 199 of DEB book.
    RAK_H = p.K_H * s.A # mol.L⁻¹ TODO: only partial root area is not involved in N uptake [@robinson1991what].
    K1_NH = half_saturation(RAK_H, SAX_H, p.K_NH) # mol.L⁻¹ Ammonia saturation
    K1_NO = half_saturation(RAK_H, SAX_H, p.K_NO) # mol.L⁻¹ Nitrate saturation
    j1_NH = half_saturation(p.X_NH, K1_NH, p.j_NH_Amax) # Arriving ammonia mols.m⁻².s⁻¹ 
    j1_NO = half_saturation(p.X_NO, K1_NO, p.j_NO_Amax) # Arriving nitrate mols.m⁻².s⁻¹

    mass_area_coef = p.w_V * p.SLA # m².mol⁻¹ = m⁻².g * g.mol⁻¹
    area = s.u.V * s.A * mass_area_coef # m² = mol * dimless * m².mol⁻¹
    J1_NH_ass = area * j1_NH # mol.s⁻¹ = m² * mols.m⁻².s⁻¹
    J_NO_ass = area * j1_NO # mol.s⁻¹ = m² * mols.m⁻².s⁻¹
    J_N_ass = J1_NH_ass + p.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end
