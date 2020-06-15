abstract type AbstractGermination end

"""
    ThresholdGermination(germination_mass)

Germination occurs past a threshhold structural mass. 

This causes a hard switch in behaviour between
"""
@columns struct ThresholdGermination{Mo} <: AbstractGermination
    # Field              | Def  | Unit | Bounds      | Log  | Description
    germination_mass::Mo | 1e-5 | mol  | (1e-10,5.0) | true | "Structural mass at germination"
end

"""
    is_germinated(formulation, o, u)

Check if germination has happened. 
The default with no formulation is that germination occurs immediately.
"""
is_germinated(o, u) = is_germinated(germination_pars(o), o, u)
is_germinated(f::Nothing, o, u) = true
is_germinated(f::ThresholdGermination, o, u) = u[:V] > f.germination_mass
