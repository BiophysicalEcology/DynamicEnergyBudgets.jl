"""
Synthesizing units bind multiple substrates to synthesize compounds, 
depending on their availability
"""
abstract type AbstractSynthesizingUnit end

"""
    synthesizing_unit(::AbstractSynthesizingUnit, v, w)

Apply a synthesizing unit formulation to substrates
`v` and `w`, returning the amount of compound.
"""
function synthesizing_unit end

"""
    ParallelComplementarySU(k)

0-parameter synthesizing unit that merges two compounds stoichiometrically.

See Ledder et al 2019. for details.

$(FIELDDOCTABLE)
"""
struct ParallelComplementarySU <: AbstractSynthesizingUnit end

synthesizing_unit(::ParallelComplementarySU, v, w) = v * w * (v + w) / (v^2 + w^2 + v * w)


"""
    MinimumRuleSU(k)

0-parameter synthesizing unit where law of the minimum controls
the production of one compound form two other compounds.

$(FIELDDOCTABLE)
"""
struct MinimumRuleSU <: AbstractSynthesizingUnit end

synthesizing_unit(::MinimumRuleSU, v, w) = min(v, w)

"""
    KfamilySU(k)

Flexible 1-parameter synthesizing unit with variable curve. Both `MinimumRuleSU`
and `ParallelComplementarySU` can be approximated with this rule.

$(FIELDDOCTABLE)
"""
@columns struct KfamilySU{K} <: AbstractSynthesizingUnit 
    k::K | 1.0 | _ | (0.0, 10.0)  | _ | "Synthesizing unit parameter. Effiency = 2^-1/k"
end

synthesizing_unit(f::KfamilySU, v, w) = (v^-f.k + w^-f.k)^(-1/f.k)

"""
    stoich_merge(Ja, Jb, ya, yb)

Merge fluxes stoichiometrically into general reserve Eab based on yield
fractions ya and yb. An unmixed proportion is returned as unmixed reserves Ea and Eb.
"""
stoich_merge(su, Jv, Jw, yv, yw) = begin
    JEvw = synthesizing_unit(su, Jv * yv, Jw * yw)
    (Jv - JEvw / yv), (Jw - JEvw / yw), JEvw
end
