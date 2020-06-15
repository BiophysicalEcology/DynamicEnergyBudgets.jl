
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
    reserve = (u[:C], u[:N])
    rel_reserve = reserve ./ u[:V]
    corrected_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r, isalive = calc_rate(su_pars(o), rel_reserve, turnover, corrected_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o), tstep(o))
    set_rate!(o, r)

    J1[:C,:ctb], J1[:N,:ctb] = non_growth_flux.(reserve, turnover, r)
    C_rej, N_rej, J1[:E,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    J[:N,:rej] = -C_rej
    J[:C,:rej] = -N_rej
    isalive
end

"""
    CatabolismCNshared(k)

1-pararameter catabolism where reserve turnover rate is the same for both reserves.
"""
@columns struct CatabolismCNshared{K} <: AbstractCatabolismCN 
    k::K  | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "Reserve turnover rate"
end

"""
    CatabolismCNshared(kC, kN)

2-pareter catabolism where reserve turnover rate can be different for C and N reserves.
"""
@columns struct CatabolismCN{KC,KN} <: AbstractCatabolismCN 
    kC::KC | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::KN | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

kC(p::CatabolismCNshared) = p.k
kN(p::CatabolismCNshared) = p.k


abstract type AbstractCatabolismCNE <: AbstractCatabolism end

catabolism!(p::AbstractCatabolismCNE, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (kC(p), kN(p), kE(p)) .* tempcorrection(v) .* shape(v)
    reserve = (u[:C], u[:N], u[:E])
    rel_reserve = reserve ./ u[:V]
    corrected_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r, isalive = calc_rate(su_pars(o), rel_reserve, turnover, corrected_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o), tstep(o))
    set_rate!(o, r)

    J1[:C,:ctb], J1[:N,:ctb], J1[:EE,:ctb] = non_growth_flux.(reserve, turnover, r)
    C_rej, N_rej, J1[:CN,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    J[:N,:rej] = -C_rej
    J[:C,:rej] = -N_rej

    J1[:E,:ctb] = J1[:EE,:ctb] + J1[:CN,:ctb] # Total catabolic flux
    set_θE!(o, J1[:EE,:ctb]/J1[:E,:ctb]) # Proportion of general reserve flux in total catabolic flux
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

"""
    growth!(o, u)

Allocates reserves to growth.
"""
growth!(o, u) = begin
    flux(o)[:V,:gro] = growth = rate(o) * u[:V]
    drain = (1/y_V_E(o)) * growth 
    product = growth_production!(o, growth)
    reserve_drain!(o, :gro, drain)
end
