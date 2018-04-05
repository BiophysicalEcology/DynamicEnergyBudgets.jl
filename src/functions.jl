# Library of common DEB functions.

export find_rate, tempcorr, save_rate!, catabolic_fluxes, stoich_merge, synthesizing_unit

"""
    area_mass_kooijman(uV, Vref, Vscaling)
Area mass scaling for plants.
(uV / Vref)^(1/3 - (-uV / Vscaling)^b) DEB book p. 134
"""
function area_mass_kooijman(uV::Float64, Vref::Float64, Vscaling::Float64)::Float64
    uV > 0.0 || return 0.0
    (uV / Vref)^(-uV / Vscaling)
end

function area_mass_linear(uV::Float64, Vref::Float64, Vscaling::Float64)::Float64
    return 1.0
end

"""
tempcorr(T, T1, A)
"""
function tempcorr(T, T1, A)
    exp(A/T1 - A/T)
end
function tempcorr(T, T1, A, L, AL)
    exp(A/T1 - A/T) * (1.0 + exp(AL/T1 - AL/L)) / (1.0 + exp(AL/T - AL/L))
end
function tempcorr(T, T1, A, L, AL, H, AH)
    exp(A/T1 - A/T) * (1.0 + exp(AL/T1 - AL/L) + exp(AH/H - AH/T1)) / 
    (1.0 + exp(AL/T - AL/L) + exp(AH/H - AH/T))
end

"""
    find_rate(t, rates::Array{Float64}, args)
Calculate rate formula.

TODO: use Roots.jl for this
"""
function find_rate(t, rates::Array{Float64}, 
                   args::Tuple{NTuple{N,Float64},NTuple{N,Float64},Vararg{Float64,V}})::Float64 where {N,V}
    local f = rate_formula
    local found = false

    # Find the largest possible rate window.
    x0, x1 = rate_window(args...) 

    # Find root with the secant method.
    y0 = f(x0, args...)
    y1 = f(x1, args...)
    i = 1
    local x::Float64
    for _ in 1:BI_MAXITER
        x = x1 - y1 * (x1 - x0) / (y1 - y0)
        if abs(x-x1) < BI_XTOL
            found = true
            break
        end
        x0 = x1
        y0 = y1
        x1 = x
        y1 = f(x1, args...)
    end

    if found
        # Save rate for plotting.
        save_rate!(t, rates, x)
        return x
    else
        error("Root not found for args: $args")
    end
end

function save_rate!(t::Float64, rates::Array{Float64}, newrate::Float64)
    t_floor = floor(Int64, t)
    if t_floor == t
        r = rates[t_floor + 1] = newrate
    end
end

"""
Get window of possible rates with CN or CNE rserves
Calculated for y_E_CH_NO and y_E_EN -> ∞
"""
function rate_window(ureserve::NTuple{2,Float64}, A_turnover::NTuple{2,Float64}, 
                     args::Vararg)::NTuple{2,Float64}
    rate_window((ureserve[1], ureserve[2], 0.0), (A_turnover[1], A_turnover[2], 0.0), args...)
end
function rate_window(ureserve::NTuple{3,Float64}, A_turnover::NTuple{3,Float64}, 
                     j_E_mai::Float64, y_E_CH_NO::Float64, y_E_EN::Float64, 
                     y_V_E::Float64, κsoma::Float64)::NTuple{2,Float64}
    # Calculate the limits of the rate_formula function when C or N → ∞
    (uEC, uEN, uE) = ureserve
    (AEC, AEN, AE) = A_turnover
    lim(y) = (uE * AE + uEC * AEC - j_E_mai/(κsoma * y))/(1/(y_V_E * κsoma * y) + uE + uEC * y)
    (lim(y_E_EN), lim(y_E_CH_NO))
end

"""
Rate formulas for E, CN or CNE reserves
"""
function rate_formula(r::Float64, ureserve::NTuple{1,Float64}, A_turnover::NTuple{1,Float64}, 
                      j_E_mai::Float64, y_V_E::Float64, κsoma::Float64)::Float64
    j_E = catabolic_fluxes(ureserve, A_turnover, r)
    return y_V_E * (κsoma * j_E - j_E_mai) - r
end
function rate_formula(r::Float64, ureserve::NTuple{2,Float64}, A_turnover::NTuple{2,Float64}, 
                      j_E_mai::Float64, y_E_CH_NO::Float64, y_E_EN::Float64, 
                      y_V_E::Float64, κsoma::Float64)::Float64
    (j_EC, j_EN) = catabolic_fluxes(ureserve, A_turnover, r)
    j_ECN = stoich_merge(j_EC * y_E_CH_NO, j_EN * y_E_EN)  
    return y_V_E * (κsoma * j_ECN - j_E_mai) - r
end
function rate_formula(r::Float64, ureserve::NTuple{3,Float64}, A_turnover::NTuple{3,Float64}, 
                      j_E_mai::Float64, y_E_CH_NO::Float64, y_E_EN::Float64, 
                      y_V_E::Float64, κsoma::Float64)::Float64
    (j_EC, j_EN, j_E) = catabolic_fluxes(ureserve, A_turnover, r)
    j_ECN = stoich_merge(j_EC * y_E_CH_NO, j_EN * y_E_EN)  
    return y_V_E * (κsoma * (j_E + j_ECN) - j_E_mai) - r
end

"""
Returns the current catabolic flux at rate r,
or the flux as a proportion of u[V], depending on ureserve values.
"""
function catabolic_fluxes(ureserve::NTuple{N,Float64}, A_turnover::NTuple{N,Float64}, 
                          r::Float64)::NTuple{N,Float64} where N
    ureserve .* (A_turnover .- r)
end

###############################################
# Stoichiometry


"""
    half_saturation(x, half, max)
Half satration curve.
"""
function half_saturation(x, half, max)
    max/(1.0 + half/x) 
end 

"""
    stoich_merge(a, b)
Merge two inputs stoichiometrically. The minimum value is limiting,
and stochasticity of pairing is simulated so that for any a, b 
stoich_merge(a, b) < a, to the limit b → ∞ where stoich_merge(a, b) = a   
"""
function stoich_merge(a::Float64, b::Float64)::Float64
    1.0/(1.0/a + 1.0/b - 1.0/(a + b))
end

"""
    synthesizing_unit(Ja, Jb, ya, yb)
Merge fluxes stoichiometrically into general reserve Eab
The remainder is returned as unmixed reserves Ea and Eb
"""
function synthesizing_unit(Ja::Float64, Jb::Float64, ya::Float64, yb::Float64
                          )::NTuple{3,Float64} 
    JEab = stoich_merge(Ja * ya, Jb * yb) 
    JEa = Ja - JEab/ya                        
    JEb = Jb - JEab/yb                     
    (JEa, JEb, JEab)
end

"""
Check if germination has happened. Independent for each structures,
although this may not make sense.
"""
function germinated(M_V, M_Vgerm::Float64)::Bool 
    M_V > M_Vgerm 
end
