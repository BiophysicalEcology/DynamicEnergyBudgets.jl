abstract type AbstractSynthesizingUnit end

"""
    ParallelComplementarySU(k)

0-parameter synthesizing unita.
"""
struct ParallelComplementarySU <: AbstractSynthesizingUnit end

"""
    MinimumRuleSU(k)

0-parameter synthesizing unit where law of the minimum controls
modilisation of multiple resources.
"""
struct MinimumRuleSU <: AbstractSynthesizingUnit end

"""
    KfamilySU(k)

1-parameter synthesizing unit
"""
@columns struct KfamilySU{K} <: AbstractSynthesizingUnit 
    k::K | 1.0 | _ | (0.0, 10.0)  | _ | "Synthesizing unit parameter. Effiency = 2^-1/k"
end

"""
    synthesizing_unit(::Type, a, b)

Merge two inputs stoichiometrically. 

See Ledder et al. TODO: add paper details once it's published.
"""
synthesizing_unit(::ParallelComplementarySU, v, w) = v * w * (v + w) / (v^2 + w^2 + v * w)
synthesizing_unit(::MinimumRuleSU, v, w) = min(v, w)
synthesizing_unit(f::KfamilySU, v, w) = (v^-f.k + w^-f.k)^(-1/f.k)

"""
    stoich_merge(Ja, Jb, ya, yb)

Merge fluxes stoichiometrically into general reserve Eab based on yeild
fractions ya and yb. An unmixed proportion is returned as unmixed reserves Ea and Eb.
"""
stoich_merge(su, Jv, Jw, yv, yw) = begin
    JEvw = synthesizing_unit(su, Jv * yv, Jw * yw)
    (Jv - JEvw / yv), (Jw - JEvw / yw), JEvw
end
