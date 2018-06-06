
function photosynthesis(f::C3Photosynthesis, o)
    o.vars.assimilation.aleaf * f.SLA * o.params.w_V
end

function photosynthesis(f::KooijmanPhotosynthesis, o)
    p = o.params
    j1_l = half_saturation(p.j_L_Amax, p.J_L_K, p.J_L_F)
    j1_c = half_saturation(p.j_C_Amax, p.K_C, p.X_C)
    j1_o = half_saturation(p.j_O_Amax, p.K_O, p.X_O)

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
    mass_area_coef = o.paramo.w_V * f.SLA
    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, f.J_L_F) * mass_area_coef
    j1_c = half_saturation(f.j_C_Amax, f.K_C, f.X_C) * mass_area_coef
    j1_o = half_saturation(f.j_O_Amax, f.K_O, f.X_O) * mass_area_coef

    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # c flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co) # dimless

    j_c_intake / (1 + bound_c + bound_o + co_l)
end

function assimilation!(o, on)
    assimilation!(o, on, o.params.assimilation, o.state)
end

function assimilation!(o, on, f::CarbonAssimilation, u::AbstractStateCNE)::Void
    p = o.params; v = o.vars
    if !germinated(u.V, p.M_Vgerm) 
        on.J[:N,:rej] = o.J[:E,:ass] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1])
        return nothing
    end

    J1_EC_ass = photosynthesis(f, o) * u.V * v.scale
    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * on.J1[:N,:rej] 
    # Merge rejected N from root and photosynthesized C into reserves
    (o.J[:C,:ass], o.J[:N,:ass], o.J[:E,:ass]) = 
        synthesizing_unit(J1_EC_ass, J1_EN_ass, p.y_E_CH_NO, p.y_E_EN)
    return nothing
end

function assimilation!(o, on, f::NH4_NO3_Assimilation, u::AbstractStateCNE)
    p = o.params
    if !germinated(u.V, p.M_Vgerm) 
        o.J[:E,:ass] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1]) 
        J[:C,:ass] = -p.y_EC_ECT * o.J1[:C,:trans]
        return nothing
    end
    (J_N_ass, J_NO_ass, J_NH_ass) = uptake_nitrogen(f, o, on)

    J1_EC_ass = -p.y_EC_ECT * on.J1[:C,:rej] # Rejected C-reserve from shoot

    θNH = J_NH_ass/J_N_ass                          # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                   # Fraction of nitrate in arriving N-flux
    y_E_CH = θNH * f.y_E_CH_NH + θNO * p.y_E_CH_NO  # Yield coefficient from C-reserve to reserve

    # Merge rejected C from shoot and uptaken N into reserves
    (o.J[:C,:ass], o.J[:N,:ass], o.J[:E,:ass]) = 
        synthesizing_unit(J1_EC_ass, J_N_ass, y_E_CH, 1/p.n_N_E)

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    o.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * o.J[:E,:ass]) * 1/p.n_N_EN
    return nothing
end
 
function uptake_nitrogen(f::Kooijman_NH4_NO3_Assimilation, o, on)
    p = o.params; v = o.vars
    K1_NH = half_saturation(f.K_NH, f.K_H * v.scale, f.X_H * on.vars.scale) # Ammonia saturation
    K1_NO = half_saturation(f.K_NO, f.K_H * v.scale, f.X_H * on.vars.scale) # Nitrate saturation
    J1_NH_ass = o.state.V * v.scale * half_saturation(f.j_NH_Amax, K1_NH, f.X_NH) # Arriving ammonia mols.mol⁻¹.s⁻¹
    J_NO_ass = o.state.V * v.scale * half_saturation(f.j_NO_Amax, K1_NO, f.X_NO) # Arriving nitrate mols.mol⁻¹.s⁻¹

    J_N_ass = J1_NH_ass + f.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.

Without a function like this you will likely be occasionally breaking the 
laws of thermodynamics by introducing negative rates.
"""
feedback!(o, f::Autophagy, u::AbstractState) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, o.vars.rate)
    autophagy = u.V * (oneunit(hs) - hs)
    o.J[:C,:gro] += autophagy/o.params.y_E_CH_NO
    o.J[:N,:gro] += autophagy/o.params.y_E_EN
    o.J[:V,:gro] -= autophagy
    nothing
end

scaling(f::KooijmanArea, uV) = begin
    uV > zero(uV) || return zero(uV)
    (uV / f.M_Vref)^(-uV / f.M_Vscaling)
end

scaling(f, uV) = uV

"""
Check if germination has happened. Independent for each organ,
although this may not make sense.
"""
germinated(M_V, M_Vgerm) = M_V > M_Vgerm 

