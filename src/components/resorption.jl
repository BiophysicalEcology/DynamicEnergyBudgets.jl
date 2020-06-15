"""
Resorption. Parameters for reabsorbtion of nutrients from structures when metabolic rates fall.
"""
abstract type AbstractResorption end


resorption!(o, u) = resorption!(resorption_pars(o), o, u)
resorption!(f::Nothing, o, u) = nothing

"""
    StructuralLossResorption(K_resorption) 

Structure is lost while reserves are retained.
"""
@columns struct StructuralLossResorption{R} <: AbstractResorption 
    K_resorption::R | 0.000001 | _    | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
end
resorption!(f::StructuralLossResorption, o, u) = begin
    o.J[:V,:res] -= resorption(f, o, u)
    nothing
end

"""
    LosslessResorption(K_resorption)

Structure is distributed back to C an N reserves without loss.
"""
@columns struct LosslessResorption{R} <: AbstractResorption 
    K_resorption::R | 0.000001 | _    | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
end
resorption!(f::AbstractResorption, o, u) = begin
    aph = resorption(f, o, u)
    o.J[:C,:res] += aph
    o.J[:N,:res] += aph * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

"""
    DissipativeResorption(r_EN_V, r_EN_V, K_resorption) 

Some structure is distributed back to C an N reserves with a proportion lost to the environment.

This model has parameters for controlling differential fractions of C an N resorption.
"""
@columns struct DissipativeResorption{R,P} <: AbstractResorption
    # Field         | Default  | Unit       | Bounds       | Log  | Description
    K_resorption::R | 0.000001 | _          | (1e-8, 1e-3) | true | "Half saturation metabolic rate for resorption of tissues."
    r_EC_V::P       | 0.0      | mol*mol^-1 | (0.0, 1.0)   | _    | "Proportion of C recovered from structure"
    r_EN_V::P       | 0.5      | mol*mol^-1 | (0.0, 1.0)   | _    | "Proportion of N recovered from structure"
end
resorption!(f::DissipativeResorption, o, u) = begin
    aph = resorption(f, o, u)
    o.J[:C,:res] += aph * f.r_EC_V
    o.J[:N,:res] += aph * f.r_EN_V * n_N_V(o)
    o.J[:V,:res] -= aph
    nothing
end

"""
    resorption(m::AbstractResorption, o::Organ, u) 

Resoption of structure `V` with shape-adjusted metabolic rate
"""
@inline resorption(m, o, u) = 
    resorption(u[:V], rate(o.vars) / shape(o), m.K_resorption)
@inline resorption(v::Number, r::Number, a::Number) = 
    v * (oneunit(r) - 1 / (oneunit(1 / r) + a / r))
