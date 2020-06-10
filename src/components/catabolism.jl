
@mix @columns struct Catabolism{MoMoD}
    # Field   | Def | Unit            | Bounds    | Log | Description
    k::MoMoD  | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "Reserve turnover rate"
end

@mix @columns struct CatabolismCN{MoMoD}
    kC::MoMoD | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::MoMoD | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

# @Catabolism struct Catabolism{} <: AbstractCatabolism end

# κEC::F      | 0.3 | _               | (0.0,1.0) | _   | "Non-processed C-reserve returned to C-reserve"
# κEN::F      | 0.3 | _               | (0.0,1.0) | _   | "Non-processed N-reserve returned to N-reserve"
 
abstract type AbstractCatabolism end

"""
    catabolism!(o, u, t::Number)

Catabolism for E, C and N, or C, N and E reserves.
Does not alter flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = catabolism!(catabolism_pars(o), o, u)

kC(p::AbstractCatabolism) = p.kC
kN(p::AbstractCatabolism) = p.kN
kE(p::AbstractCatabolism) = p.k

abstract type AbstractCatabolismCN <: AbstractCatabolism end

catabolism!(p::AbstractCatabolismCN, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (kC(p), kN(p)) .* tempcorrection(v) .* shape(v)
    reserve = (u.C, u.N)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r, found = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    if !found
        println("Root for rate not found at t - $(tstep(o))")
        return false # dead
    elseif r < zero(r) 
        println("Rate is less than zero at t - $(tstep(o))")
        return false # dead
    else
        set_rate!(o, r)
    end

    J1[:C,:ctb], J1[:N,:ctb] = reserve .* (turnover .- r)
    J1[:C,:rej], J1[:N,:rej], J1[:E,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    true
end

"""
    CatabolismCNshared(k)

1-pararameter catabolism where reserve turnover rate is the same for both reserves.
"""
@Catabolism struct CatabolismCNshared{} <: AbstractCatabolismCN end

"""
    CatabolismCNshared(kC, kN)

2-pareter catabolism where reserve turnover rate can be different for C and N reserves.
"""
@CatabolismCN struct CatabolismCN{} <: AbstractCatabolismCN end

kC(p::CatabolismCNshared) = p.k
kN(p::CatabolismCNshared) = p.k


abstract type AbstractCatabolismCNE <: AbstractCatabolism end

catabolism!(p::AbstractCatabolismCNE, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (kC(p), kN(p), kE(p)) .* tempcorrection(v) .* shape(v)
    reserve = (u.C, u.N, u.E)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    set_rate!(o, r)
    r < zero(r) && return false

    J1[:C,:ctb], J1[:N,:ctb], J1[:EE,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:CN,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))

    J1[:E,:ctb] = J1[:EE,:ctb] + J1[:CN,:ctb] # Total catabolic flux
    set_θE!(o, J1[:EE,:ctb]/J1[:E,:ctb]) # Proportion of general reserve flux in total catabolic flux
    true
end

"""
    CatabolismCNshared(k, kC, kN)

3-pareter catabolism where reserve turnover rate can be different for C, N and E reserves.
"""
@Catabolism @CatabolismCN struct CatabolismCNE{} <: AbstractCatabolismCNE end


abstract type AbstractRate end

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
    (j_E,) = rel_reserve .* (turnover .- r)
    y_V_E * (κsoma * j_E - j_E_mai) - r

    j_E = reserve * (turnover - y_V_E * (κsoma * j_E - j_E_mai))
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
    r = y_V_E * (κsoma * j_E - j_E_mai)
end

"""
    growth!(o, u)

Allocates reserves to growth.
"""
growth!(o, u) = begin
    flux(o)[:V,:gro] = growth = rate(o) * u.V
    drain = (1/y_V_E(o)) * growth 
    product = growth_production!(o, growth)
    reserve_drain!(o, Val(:gro), drain)
    # loss = drain - growth - product
    # reserve_loss!(o, loss)
    # conversion_loss!(o, growth, n_N_V(o))
end
