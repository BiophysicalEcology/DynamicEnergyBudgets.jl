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

$(FIELDDOCTABLE)
"""
@columns struct CatabolismCN{KC,KN} <: AbstractCatabolismCN 
    kC::KC | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::KN | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

"""
    CatabolismCNshared(k)

1-pararameter catabolism where reserve turnover rate is the same for both reserves.

$(FIELDDOCTABLE)
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

$(FIELDDOCTABLE)
"""
@columns struct CatabolismCNE{KE,KC,KN} <: AbstractCatabolismCNE
    k::KE  | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "Reserve turnover rate"
    kC::KC | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "C-reserve turnover rate"
    kN::KN | 0.2 | mol*mol^-1*d^-1 | (0.0,1.0) | _   | "N-reserve turnover rate"
end

