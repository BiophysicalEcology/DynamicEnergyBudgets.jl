" Resorption parameters. These modfy state based on state. "
abstract type AbstractResorption end

@mix @columns struct MixinResorption{Mo}
    # Field         | Default  | Unit | Pror           | Limits       | Log  | Description
    K_autophagy::Mo | 0.000001 | _    | Beta(2.0, 2.0) | [1e-8, 1e-3] | true | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end
" esorption. Parameters for self reabsorbtion when metabolic rates fall "
@MixinResorption struct LosslessResorption{} <: AbstractResorption end
@MixinResorption struct StructuralLossResorption{} <: AbstractResorption end

@MixinResorption struct DissipativeResorption{MoMo} <: AbstractResorption
    r_EC_V::MoMo    | 0.0      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _    | "Recovery of C from structure"
    r_EN_V::MoMo    | 0.5      | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _    | "Recovery of N from structure"
end

"""
Function to apply resorption on growth the process, such as autopagy in resource shortage.
"""
resorption!(o, u) = resorption!(resorption_pars(o), o, u)
resorption!(f::Nothing, o, u) = nothing
resorption!(f::LosslessResorption, o, u) = begin
    o.J[:V,:fbk] -= autophagy(f, o, u)
    nothing
end
resorption!(f::StructuralLossResorption, o, u) = begin
    aph = autophagy(f, o, u)
    o.J[:C,:fbk] += aph
    o.J[:N,:fbk] += aph * n_N_V(o)
    o.J[:V,:fbk] -= aph
    nothing
end
resorption!(f::DissipativeResorption, o, u) = begin
    aph = autophagy(f, o, u)
    o.J[:C,:fbk] += aph * f.r_EC_V
    o.J[:N,:fbk] += aph * f.r_EN_V * n_N_V(o)
    o.J[:V,:fbk] -= aph
    nothing
end

autophagy(f, o, u) = autophagy(u.V, rate(o.vars), f.K_autophagy)
autophagy(v::Number, r::Number, a::Number) = v * (oneunit(r) - 1 / (oneunit(1/r) + a / r))
