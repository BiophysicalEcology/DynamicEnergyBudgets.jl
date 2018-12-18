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

tempcorr(o::Organ) = tempcorr(temp(o.vars), tempcorr_pars(o))
tempcorr(t, tc::Nothing) = 1.0
tempcorr(t, tc::TempCorr) = tempcorr(t |> K, tc.reftemp, tc.arrtemp)
tempcorr(t, tc::TempCorrLower) =
    tempcorr(t |> K, tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower)
tempcorr(t, tc::TempCorrLowerUpper) =
    tempcorr(t |> K, tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower, tc.upperbound, tc.arrupper)

"""
Calculate rate formula.
"""
function calc_rate(::FZeroRate, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)
    # Find the type of `rate` so that diff and units work with atol etc
    one_r = oneunit(eltype(turnover))
    # Get rate with a zero finder
    let m=m, turnover=turnover, j_E_mai=j_E_mai, y_E_Ea=y_E_Ea, y_E_Eb=y_E_Eb, y_V_E=y_V_E, κsoma=κsoma
        find_zero((-2one_r, 1one_r), Secant(), one_r*1e-10, 100) do r
            rate_formula(r, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)
        end
    end
end

"""
    rate_formula(r, ureserve::NTuple, turnover::NTuple, j_E_mai, y_V_E, κsoma)
Rate formulas for E, CN or CNE reserves
"""
rate_formula(r, su, rel_reserve::NTuple{1}, turnover::NTuple{1},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    (j_E,) = non_growth_flux.(rel_reserve, turnover, r)
    y_V_E * (κsoma * j_E - j_E_mai) - r

    j_E = reserve * (turnover - y_V_E * (κsoma * j_E - j_E_mai))
end
rate_formula(r, su, rel_reserve::NTuple{2}, turnover::NTuple{2},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb = non_growth_flux.(rel_reserve, turnover, r)
    j_E = synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    κsoma * j_E - j_E_mai - r
end
rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = non_growth_flux.(rel_reserve, turnover, r)
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    κsoma * j_E - j_E_mai - r
end


calc_rate(::SimpleRate, su, rel_reserve::NTuple{2}, turnover::NTuple{2}, 
          j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb = rel_reserve .* turnover
    j_E = synthesizing_unit(j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    j_E - j_E_mai
end
calc_rate(::SimpleRate, su, rel_reserve::NTuple{3}, turnover::NTuple{3}, 
          j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = rel_reserve .* turnover
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    j_E - j_E_mai
end


"""
    non_growth_flux(ureserve, turnover, r)
Returns the current non_growth_flux flux at rate r,
or the flux as a proportion of u[V], depending on ureserve values.
"""
non_growth_flux(reserve, turnover, r) = reserve * (turnover - r)

"""
    half_saturation(max, half, x)
Half satration curve.
"""
half_saturation(max, half, x) = max/(oneunit(half/x) + half/x)

"""
    synthesizing_unit(::Type, a, b)
Merge two inputs stoichiometrically. The minimum value is limiting,
and stochasticity of pairing is simulated so that for any a, b
stoich_merge(a, b) < a, to the limit b → ∞ where stoich_merge(a, b) = a

SU(v, w) = 1/(1/v + 1/w - 1/(v + w))
"""
synthesizing_unit(o::Organ, v, w) = synthesizing_unit(su_pars(o), v, w) 
synthesizing_unit(::ParallelComplementarySU, v, w) = v * w * (v + w) / (v^2 + w^2 + v * w)
synthesizing_unit(::MinimumRuleSU, v, w) = min(v, w)
synthesizing_unit(f::KfamilySU, v, w) = 1/(v^-f.k + w^-f.k)^(1/f.k)
