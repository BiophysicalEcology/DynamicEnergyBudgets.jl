
@mix @columns struct Catabolism{MoMoD}
    k::MoMoD  | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "Reserve turnover rate"
end

@mix @columns struct CatabolismCN{MoMoD}
    kC::MoMoD | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "C-reserve turnover rate"
    kN::MoMoD | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "N-reserve turnover rate"
end

abstract type AbstractCatabolism end
abstract type AbstractCatabolismCN <: AbstractCatabolism end
abstract type AbstractCatabolismCNE <: AbstractCatabolism end
@Catabolism struct Catabolism{} <: AbstractCatabolism end
@Catabolism struct CatabolismCNshared{} <: AbstractCatabolismCN end
@CatabolismCN struct CatabolismCN{} <: AbstractCatabolismCN end
@Catabolism @CatabolismCN struct CatabolismCNE{} <: AbstractCatabolismCNE end

# κEC::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]  | _ | "Non-processed C-reserve returned to C-reserve"
# κEN::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]  | _ | "Non-processed N-reserve returned to N-reserve"

kC(p::CatabolismCNshared) = p.k
kN(p::CatabolismCNshared) = p.k
kC(p::AbstractCatabolism) = p.kC
kN(p::AbstractCatabolism) = p.kN
kE(p::AbstractCatabolism) = p.k

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not alter flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = catabolism!(catabolism_pars(o), o, u)
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

    J1[:C,:ctb], J1[:N,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:E,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    true
end
