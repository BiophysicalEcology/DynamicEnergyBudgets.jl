abstract type AbstractSU end

struct ParallelComplementarySU <: AbstractSU end

struct MinimumRuleSU <: AbstractSU end

@columns struct KfamilySU{K} <: AbstractSU 
    k::K | 1.0 | _ | Gamma(2.0, 2.0) | [0.0, 10.0]  | _ | "Synthesizing unit parameter. Effiency = 2^-1/k"
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
stoich_merge(su, Ja, Jb, ya, yb) = begin
    JEab = synthesizing_unit(su, Ja * ya, Jb * yb)
    Ja1 = Ja - JEab/ya
    Jb1 = Jb - JEab/yb
    (Ja1, Jb1, JEab)
end
