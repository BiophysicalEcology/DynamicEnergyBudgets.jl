" Temperature correction parameters"
abstract type AbstractTempCorr{T} end

# Temperature mixins
@mix @columns struct Tbase{T}
    # Field       | Default | Unit | Prior               | Limits              | Description
    reftemp::T    | 310.0   | K | Gamma(310.0, 1.0)   | [273.0, 325.0]     | _ | "Reference temperature for all rate parameters"
    arrtemp::T    | 2000.0  | K | Gamma(2000.0, 1.0)  | [200.0, 4000.0]    | _ | "Arrhenius temperature"
end

@mix @columns struct Tlow{T}
    lowerbound::T | 280.0   | K | Gamma(280.0, 1.0)   | [273.0, 325.0]     | _ | "Lower boundary of tolerance range"
    arrlower::T   | 20000.0 | K | Gamma(20000.0, 1.0) | [2000.0, 40000.0]  | _ | "Arrhenius temperature for lower boundary"
end

@mix @columns struct Tup{T}
    upperbound::T | 315.0   | K | Gamma(315.0, 1.0)   | [273.0, 325.0]     | _ | "Upper boundary of tolerance range"
    arrupper::T   | 70000.0 | K | Gamma(70000, 1.0)   | [7000.0, 140000.0] | _ | "Arrhenius temperature for upper boundary"
end

# Temperature types
" Simple temperature correction parameters "
@Tbase struct TempCorr{T} <: AbstractTempCorr{T} end
" Temperature correction with lower boudn parameters"
@Tbase @Tlow struct TempCorrLower{T} <: AbstractTempCorr{T} end
" Temperature correction with lower and upper bound parameters"
@Tbase @Tlow @Tup struct TempCorrLowerUpper{T} <: AbstractTempCorr{T} end


tempcorr(t, tc::Nothing) = 1.0
tempcorr(t, tc::TempCorr) = tempcorr(t |> K, tc.reftemp, tc.arrtemp)
tempcorr(t, tc::TempCorrLower) =
    tempcorr(t |> K, tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower)
tempcorr(t, tc::TempCorrLowerUpper) =
    tempcorr(t |> K, tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower, tc.upperbound, tc.arrupper)

"""
    tempcorr(T, T1, A, [L, AL,] [H, AH])
DEB tempcorr function. Uses lower and uppper bounds if they are supplied.
Temperatures all in Kelvins.
"""
tempcorr(t, tref, a) = exp(a/tref - a/t)
tempcorr(t, tref, a, l, al) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l)) / (1.0 + exp(al/t - al/l))
tempcorr(t, tref, a, l, al, h, ah) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l) + exp(ah/h - ah/tref)) /
    (1.0 + exp(al/t - al/l) + exp(ah/h - ah/t))
