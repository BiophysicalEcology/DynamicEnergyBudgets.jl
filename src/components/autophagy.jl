" State feedback parameters. These modfy state based on state. "
abstract type AbstractStateFeedback end

@mix @columns struct KAutophagy{Mo}
    # Field         | Default  | Unit       | Pror           | Limits       | Log  | Description
    K_autophagy::Mo | 0.000001 | mol        | Beta(2.0, 2.0) | [1e-8, 1e-3] | true | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end
" Autophagy. Parameters for self reabsorbtion when metabolic rates fall "
@KAutophagy struct LosslessAutophagy{Mo} <: AbstractStateFeedback end
@KAutophagy struct StructuralLossAutophagy{Mo} <: AbstractStateFeedback end

@KAutophagy struct DissipativeAutophagy{Mo,MoMo} <: AbstractStateFeedback
    r_EC_V::MoMo    | 0.0      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _    | "Recovery of C from structure"
    r_EN_V::MoMo    | 0.5      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _    | "Recovery of N from structure"
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.
"""
feedback!(o, u) = feedback!(feedback_pars(o), o, u)
feedback!(f::Nothing, o, u) = nothing
feedback!(f, o, u) = begin
    c, n, v = feedback(f, o, u)
    o.J[:C,:fbk] += c
    o.J[:N,:fbk] += n
    o.J[:V,:fbk] += v
    nothing
end
feedback(f::LosslessAutophagy, o, u) = begin
    aph = autophagy(f, o, u)
    c, n, v = aph, aph * n_N_V(o), -aph
end
feedback(f::StructuralLossAutophagy, o, u) = begin
    aph = autophagy(f, o, u)
    c, n, v = zero(aph), zero(aph), -aph
end
feedback(f::DissipativeAutophagy, o, u) = begin
    aph = autophagy(f, o, u)
    c, n, v = aph * f.r_EC_V, aph * f.r_EN_V * n_N_V(o), -aph
end

autophagy(f, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    u.V * (oneunit(hs) - hs)
end
