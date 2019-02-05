" State feedback parameters. These modfy state based on state. "
abstract type AbstractStateFeedback end

@mix @columns struct KAutophagy{Mo}
    K_autophagy::Mo | 0.000001 | mol  | Beta(2.0, 2.0)  | [0.0000001, 0.0001] | _ | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end
" Autophagy. Parameters for self reabsorbtion when metabolic rates fall "
@KAutophagy struct LosslessAutophagy{Mo} <: AbstractStateFeedback end

@KAutophagy struct DissipativeAutophagy{Mo,MoMo} <: AbstractStateFeedback
    r_EC_V::MoMo    | 0.2      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0] | _ | "Recovery of C from structure"
    r_EN_V::MoMo    | 0.6      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0] | _ | "Recovery of N from structure"
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.
"""
feedback!(o, u) = feedback!(feedback_pars(o), o, u)
feedback!(f::Nothing, o, u) = nothing
feedback!(f::LosslessAutophagy, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    o.J[:C,:fbk] += aph
    o.J[:N,:fbk] += aph * n_N_V(o)
    o.J[:V,:fbk] -= aph
    nothing
end
feedback!(f::DissipativeAutophagy, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    aC = aph * f.r_EC_V
    aN = aph * f.r_EN_V * n_N_V(o)
    o.J[:C,:fbk] += aC
    o.J[:N,:fbk] += aN
    o.J[:V,:fbk] -= aph
    nothing
end
