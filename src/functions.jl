# Library of common DEB functions.

"""
    tempcorr(T, T1, A, [L, AL,] [H, AH])
DEB tempcorr function. Uses lower and uppper bounds if they are supplied.
Temperatures all in Kelvins.
"""
tempcorr(t, t1, a) = exp(a/t1 - a/t)
tempcorr(t, t1, a, l, al) = 
    exp(a/t1 - a/t) * (1.0 + exp(al/t1 - al/l)) / (1.0 + exp(al/t - al/l))
tempcorr(t, t1, a, l, al, h, ah) = 
    exp(a/t1 - a/t) * (1.0 + exp(al/t1 - al/l) + exp(ah/h - ah/t1)) / 
    (1.0 + exp(al/t - al/l) + exp(ah/h - ah/t))
tempcorr(t, tc::Void) = 1.0
tempcorr(t, tc::TempCorr) = tempcorr(t |> u"K", tc.reftemp, tc.arrtemp)
tempcorr(t, tc::TempCorrLower) = tempcorr(t |> u"K", tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower)
tempcorr(t, tc::TempCorrLowerUpper) = tempcorr(t |> u"K", tc.reftemp, tc.arrtemp, tc.lowerbound, tc.arrlower, tc.upperbound, tc.arrupper)

"""
    rate_formula(r, ureserve::NTuple, A_turnover::NTuple, j_E_mai, y_V_E, κsoma)
Rate formulas for E, CN or CNE reserves
"""
function rate_formula(r, ureserve::NTuple{3}, A_turnover::NTuple{3}, 
                      j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)
    (j_EC, j_EN, j_E) = catabolic_flux(ureserve, A_turnover, r)
    j_ECN = synthesizing_unit(j_EC * y_E_CH_NO, j_EN * y_E_EN)  
    y_V_E * (κsoma * (j_E + j_ECN) - j_E_mai) - r
end

function rate_bracket(ureserve::NTuple{3}, A_turnover::NTuple{3}, 
                     j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)::NTuple{2}
    # Calculate the limits of the rate_formula function when C or N → ∞
    (uEC, uEN, uE) = ureserve
    (AEC, AEN, AE) = A_turnover
    lim(y) = (uE * AE + uEC * AEC - j_E_mai/(κsoma * y))/(1/(y_V_E * κsoma * y) + uE + uEC * y)
    (lim(y_E_EN), lim(y_E_CH_NO))
end

"""
    catabolic_flux(ureserve, A_turnover, r)
Returns the current catabolic flux at rate r,
or the flux as a proportion of u[V], depending on ureserve values.
"""
catabolic_flux(ureserve, A_turnover, r) = ureserve .* (A_turnover .- r)

"""
    half_saturation(max, half, x)
Half satration curve.
"""
half_saturation(max, half, x) = max/(oneunit(half/x) + half/x) 

"""
    synthesizing_unit(a, b)
Merge two inputs stoichiometrically. The minimum value is limiting,
and stochasticity of pairing is simulated so that for any a, b 
stoich_merge(a, b) < a, to the limit b → ∞ where stoich_merge(a, b) = a   
"""
synthesizing_unit(a, b) = 1.0/(1.0/a + 1.0/b - 1.0/(a + b))
