abstract type AbstractRate end

struct SimpleRate <: AbstractRate end

struct FZeroRate <: AbstractRate end

"""
Calculate rate formula.
"""
function calc_rate(::FZeroRate, su, rel_reserve, turnover, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma)
    # Find the type of `rate` so that diff and units work with atol etc
    one_r = oneunit(eltype(turnover))
    # Get rate with a zero finder
    let m=m, turnover=turnover, j_E_mai=j_E_mai, y_E_Ea=y_E_Ea, y_E_Eb=y_E_Eb, y_V_E=y_V_E, κsoma=κsoma
        find_zero((-2one_r, 1one_r), Secant(), one_r*1e-10, 200) do r
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
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = non_growth_flux.(rel_reserve, turnover, r)
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    r = y_V_E * (κsoma * j_E - j_E_mai)
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

