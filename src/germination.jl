abstract type AbstractGermination end

@columns struct Germination{Mo} <: AbstractGermination
    M_Vgerm::Mo | 0.0 | mol             | Gamma(2.0, 2.0) | [0.0,1.0]    | _ | "Structural mass at germination"
end

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
is_germinated(o, u) = is_germinated(germination_pars(o), o, u)
is_germinated(g::Nothing, o, u) = true
is_germinated(g::Germination, o, u) = u.V > g.M_Vgerm
