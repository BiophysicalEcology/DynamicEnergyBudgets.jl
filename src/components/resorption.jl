"""
Resorption. Parameters for reabsorbtion of nutrients from structures when metabolic rates fall.
"""
abstract type AbstractResorption end

@mix @columns struct MixinResorption{Mo}
    # Field          | Default  | Unit | Bounds       | Log  | Description
    K_resorption::Mo | 0.000001 | _    | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
end

resorption!(o, u) = resorption!(resorption_pars(o), o, u)
resorption!(f::Nothing, o, u) = nothing

"""
    LosslessResorption(K_resorption) 

Complete nutrient resorption without losses to the environment.

# TODO this seems erroneous 
"""
@MixinResorption struct LosslessResorption{} <: AbstractResorption end
resorption!(f::LosslessResorption, o, u) = begin
    o.J[:V,:res] -= resorption(f, o, u)
    nothing
end

"""
    StructuralLossResorption(K_resorption)

Structure is distributed back to C an N reserves without loss.
"""
@MixinResorption struct StructuralLossResorption{} <: AbstractResorption end
resorption!(f::StructuralLossResorption, o, u) = begin
    aph = resorption(f, o, u)
    o.J[:C,:res] += aph
    o.J[:N,:res] += aph * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

"""
    DissipativeResorption(r_EN_V, r_EN_V, K_resorption) 

Some structure is distributed back to C an N reserves with a proportion lost to the environment.
"""
@MixinResorption struct DissipativeResorption{MoMo} <: AbstractResorption
    # Field         | Default  | Unit       | Bounds      | Log  | Description
    r_EC_V::MoMo    | 0.0      | mol*mol^-1 | (0.0, 1.0)  | _    | "Proportion of C recovered from structure"
    r_EN_V::MoMo    | 0.5      | mol*mol^-1 | (0.0, 1.0)  | _    | "Proportion of N recovered from structure"
end
resorption!(f::DissipativeResorption, o, u) = begin
    aph = resorption(f, o, u)
    o.J[:C,:res] += aph * f.r_EC_V
    o.J[:N,:res] += aph * f.r_EN_V * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

@inline resorption(f, o, u) = resorption(u.V, rate(o.vars) / shape(o), f.K_resorption)
@inline resorption(v::Number, r::Number, a::Number) = v * (oneunit(r) - 1 / (oneunit(1/r) + a / r))
