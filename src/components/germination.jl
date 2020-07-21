"""
Abstract supertype for models that determine the time of germination
from current state variable.
"""
abstract type AbstractGermination end

"""
    ThresholdGermination(germination_mass)

Germination occurs past a threshhold structural mass. 

This causes a hard switch in behaviour between

$(FIELDDOCTABLE)
"""
@columns struct ThresholdGermination{Mo} <: AbstractGermination
    # Field              | Def  | Unit | Bounds      | Log  | Description
    germination_mass::Mo | 1e-5 | mol  | (1e-10,5.0) | true | "Structural mass at germination"
end

"""
    isgerminated(formulation, o, u)

Check if germination has happened. 
The default with no formulation is that germination occurs immediately.
"""
isgerminated(o, u) = isgerminated(germination_pars(o), o, u)
isgerminated(f::Nothing, o, u) = true
isgerminated(f::ThresholdGermination, o, u) = u[:V] > f.germination_mass
