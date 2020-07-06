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
    y_V_E * (κsoma * j_E - j_E_mai) - r
end
rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3},
             j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
    j_Ea, j_Eb, j_E = rel_reserve .* (turnover .- r)
    j_E += synthesizing_unit(su, j_Ea * y_E_Ea, j_Eb * y_E_Eb)
    y_V_E * (κsoma * j_E - j_E_mai) - r
end


abstract type AbstractCatabolism end

"""
    catabolism!(o, u)

Where `o` is an `Organ` and `u` its state variables

Catabolism for E, C and N, or C, N and E reserves.
"""
catabolism!(o, u) = catabolism!(catabolism_pars(o), o, u)

kC(p::AbstractCatabolism) = p.kC
kN(p::AbstractCatabolism) = p.kN
kE(p::AbstractCatabolism) = p.k

abstract type AbstractCatabolismCN <: AbstractCatabolism end

catabolism!(p::AbstractCatabolismCN, o, u) = begin
    v, J = vars(o), flux(o)
    turnover = (kC(p), kN(p)) .* tempcorrection(v) .* scaling(v)
    rel_reserve = (u[:C], u[:N]) ./ u[:V]
    corrected_j_E_mai = j_E_mai(o) * tempcorrection(v) 

    # Calculate the growth rate
    r, alive = calc_rate(su_pars(o), rel_reserve, turnover, corrected_j_E_mai, 
                         y_E_C(o), y_E_N(o), y_V_E(o), κsoma(o), tstep(o))
    set_rate!(o, r)

    # Catabolise reserve depending on reserve state and turnover minus the growth rate
    C_ctb, N_ctb = (u[:C], u[:N]) .* (turnover .- r)
    C_rej, N_rej, E_ctb = stoich_merge(su_pars(o), C_ctb, N_ctb, y_E_C(o), y_E_N(o))
    # Set catabolism var for use in other components
    set_E_ctb!(o, E_ctb)
    # Set rejected catabolised reserves to be translocated, or lost
    J[:C, :rej] = -C_rej
    J[:N, :rej] = -N_rej
    alive
end

"""
    CatabolismCN(kC, kN)

2-pareter catabolism where reserve turnover rate can be different for C and N reserves.
"""
@columns struct CatabolismCN{KC,KN} <: AbstractCatabolismCN 
    kC::KC | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::KN | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

"""
    CatabolismCNshared(k)

1-pararameter catabolism where reserve turnover rate is the same for both reserves.
"""
@columns struct CatabolismCNshared{K} <: AbstractCatabolismCN 
    k::K  | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "Reserve turnover rate"
end

kC(p::CatabolismCNshared) = p.k
kN(p::CatabolismCNshared) = p.k


abstract type AbstractCatabolismCNE <: AbstractCatabolism end

catabolism!(p::AbstractCatabolismCNE, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (kC(p), kN(p), kE(p)) .* tempcorrection(v) .* scaling(v)
    reserve = (u[:C], u[:N], u[:E])
    rel_reserve = reserve ./ u[:V]
    corrected_j_E_mai = j_E_mai(o) * tempcorrection(o) 

    r, isalive = calc_rate(su_pars(o), rel_reserve, turnover, corrected_j_E_mai(o), y_E_C(o), y_E_N(o), y_V_E(o), κsoma(o), tstep(o))
    set_rate!(o, r)

    J1_C_ctb, J1_N_ctb, J_EE_ctb = reserve .* (turnover .- r)
    C_rej, N_rej, J1_CN_ctb = stoich_merge(su_pars(o), J_C_ctb, J1_N_ctb, y_E_C(o), y_E_N(o))
    E_ctb = J1_EE_ctb + J1_CN_ctb # Total catabolic flux

    # Set rejected flux
    J[:N,:rej] = -C_rej
    J[:C,:rej] = -N_rej
    # Set catabolism vars for use in other components
    set_θE!(o, J_EE_ctb / J_E_ctb) # Proportion of general reserve flux in total catabolic flux
    set_E_ctb!(o, J_E_ctb)
    isalive
end

"""
    CatabolismCNshared(k, kC, kN)

3-pareter catabolism where reserve turnover rate can be different for C, N and E reserves.
"""
@columns struct CatabolismCNE{KE,KC,KN} <: AbstractCatabolismCNE
    k::KE  | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "Reserve turnover rate"
    kC::KC | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::KN | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

