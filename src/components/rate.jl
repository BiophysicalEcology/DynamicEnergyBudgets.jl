
"""
Calculate rate formula.
"""
function calc_rate(su, rel_reserve::Tuple, turnover::Tuple, 
                   j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma, tstep)
    # Find the type of `rate` so that diff and units work with atol etc
    one_r = oneunit(eltype(turnover))
    # Get rate with a zero finder
    r, state = let turnover=turnover, j_E_mai=j_E_mai, y_E_Ea=y_E_Ea, y_E_Eb=y_E_Eb, y_V_E=y_V_E, κsoma=κsoma
        atol = one_r * 1e-10; maxiters = 200
        findzero(Secant(-2one_r, 1one_r), atol, maxiters) do r
            rate_formula(r, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)
        end
    end

    if state == DEAD
        @warn "Root for rate not found at t = $tstep"
        return zero(r), DEAD
    elseif r < zero(r) 
        @warn "Rate is less than zero at t = $tstep"
        return zero(r), DEAD
    else 
        return r, ALIVE
    end
end

"""
    rate_formula(r, ureserve::NTuple, turnover::NTuple, j_E_mai, y_V_E, κsoma)

Rate formulas for E, CN or CNE reserves
"""
rate_formula(r, su, rel_reserve::NTuple{1}, turnover::NTuple{1},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    (j_E,) = rel_reserve .* (turnover .- r)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
rate_formula(r, su, rel_reserve::NTuple{2}, turnover::NTuple{2},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb = rel_reserve .* (turnover .- r)
    j_E = synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    r = y_V_E * (κsoma * j_E - j_E_mai)
end
rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = rel_reserve .* (turnover .- r)
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
