abstract type AbstractGermination end

@columns struct ThresholdGermination{Mo} <: AbstractGermination
    # Field              | Def  | Unit | Pror            | Limits      | Log  | Description
    germination_mass::Mo | 1e-5 | mol  | Gamma(2.0, 2.0) | [1e-10,5.0] | true | "Structural mass at germination"
end

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
is_germinated(o, u) = is_germinated(germination_pars(o), o, u)
is_germinated(f::Nothing, o, u) = true
is_germinated(f::ThresholdGermination, o, u) = u.V > g.germination_mass
