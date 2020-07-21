"""
Resorption. Parameters for reabsorbtion of nutrients from structures when metabolic rates fall.
"""
abstract type AbstractResorption end

K_resorption(p::AbstractResorption) = p.K_resorption

resorption!(o, u) = resorption!(resorption_pars(o), o, u)
resorption!(p::Nothing, o, u) = nothing

"""
    StructuralLossResorption(K_resorption) 

Structure is lost while reserves are retained.

$(FIELDDOCTABLE)
"""
@columns struct StructuralLossResorption{R} <: AbstractResorption 
    K_resorption::R | 0.000001 | _    | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
end

resorption!(p::StructuralLossResorption, o, u) = begin
    o.J[:V,:res] -= resorption(p, o, u)
    nothing
end

"""
    LosslessResorption(K_resorption)

Structure is distributed back to C an N reserves without loss.

$(FIELDDOCTABLE)
"""
@columns struct LosslessResorption{R} <: AbstractResorption 
    K_resorption::R | 0.000001 | _    | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
end

resorption!(p::LosslessResorption, o, u) = begin
    aph = resorption(p, o, u)
    o.J[:C,:res] += aph
    o.J[:N,:res] += aph * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

"""
    DissipativeResorption(r_EN_V, r_EN_V, K_resorption) 

Some structure is distributed back to C an N reserves with a proportion lost to the environment.

This model has parameters for controlling differential fractions of C an N resorption.

$(FIELDDOCTABLE)
"""
@columns struct DissipativeResorption{R,P} <: AbstractResorption
    # Field         | Default  | Unit       | Bounds       | Log  | Description
    K_resorption::R | 0.000001 | _          | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
    r_EC_V::P       | 0.0      | mol*mol^-1 | (0.0, 1.0)   | _    | "Proportion of C recovered from structure"
    r_EN_V::P       | 0.5      | mol*mol^-1 | (0.0, 1.0)   | _    | "Proportion of N recovered from structure"
end

resorption!(p::DissipativeResorption, o, u) = begin
    aph = resorption(m, o, u)
    o.J[:C,:res] += aph * p.r_EC_V
    o.J[:N,:res] += aph * p.r_EN_V * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

"""
    resorption(p::AbstractResorption, o::Organ, u) 

Resoption of structure `V` with shape-adjusted metabolic rate
"""
@inline resorption(p::AbstractResorption, o, u) = 
    resorption(u[:V], rate(o) / scaling(o), K_resorption(p))

@inline resorption(v::Number, r::Number, a::Number) = 
    v * (oneunit(r) - 1 / (oneunit(1 / r) + a / r))
