import MechanisticModels.runmodel!

"""
    runmodel!(du, settings, u, t)
DEBSettings method for MechanisticModels.jl api.
Applies environment and runs the DEB model.
"""
function runmodel!(du, settings::S, u, t::Number)::Void where S<:DEBSettings
    ss = settings.structures
    split_state!(ss, 0, u)

    if settings.use_environment
        settings.apply_environment!(settings, t)
    end

    # Store intermediate things for plotting. Slower and more memory intensive.
    if settings.save_intermediate
        set_current_flux!, ss, floor(Int, t)
    end

    deb_model!(settings, t)
    sum_flux!(du, settings)
    return nothing
end


"""
    deb_model!(settings, t)
A generalised multi-reserve, multi-structure Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N structures.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
function deb_model!(settings::S, t::Number)::Void where S
    structures = settings.structures
    swapped = (Base.tail(structures)..., structures[1])

    apply(metabolism!, structures, t)
    if length(structures) > 1
        # FIXME what about rejected reserves in 1 structure organism?
        apply(translocation!, structures, swapped)
    end
    apply(assimilation!, structures, swapped, settings)

    return nothing
end

"""
Metabolism is an identical process for all structures, with potentially
different parameters or area and rate functions.
"""
function metabolism!(s::S, t::Number)::Void where S<:DEBStructure
    s.A = s.functions.area(s.u.V, s.params.M_Vref, s.params.M_Vscaling)
    (θE, r) = catabolism!(s, s.u, t)
    dissipation!(s, s.u, θE, r)
    # autophagy!(u[i], params, J[i][i], r)
    return nothing
end


"""
Some rejected reserve is translocated.
But this grouping is somewhat unsatisfying.
"""
function translocation!(s::S1, sn::S2)::Void where {S1,S2}
    reuse_rejected_reserve!(s, sn)
    translocate!(s, sn, s.u)
    return nothing
end


"""
Run a structure-specific assimilation function.
"""
function assimilation!(s::S1, sn::S2, p::P)::Void where {S1,S2,P}
    s.functions.assim(s, sn, s.u, p)
    return nothing
end


"""
Catabolism for E, C and N, or C, N and E reserves.
"""
function catabolism!(s::S, u::AbstractState, t::Float64)::NTuple{2,Float64} where S
    return (0.0, 0.0)
end
function catabolism!(s::S, u::AbstractStateE, t)::NTuple{2,Float64} where S
    p = s.params
    J1 = s.J1
    Aturnover = p.k_E * s.A
    m = u.E / u.V
    r = rate(t, s.rates, (m, Aturnover, p.j_E_mai, p.y_V_E, p.κsoma))
    J1[:E,:cat] = catabolic_fluxes(ureserve, Aturnover, r)
    return (0.0, r)
end
function catabolism!(s::S, u::AbstractStateCN, t::Float64)::NTuple{2,Float64} where S
    p = s.params
    J1 = s.J1
    Aturnover = (p.k_EC, p.k_EN) .* s.A
    ureserve = (u.C, u.N)
    m = ureserve ./ u.V
    r = find_rate(t, s.rates, (m, Aturnover, p.j_E_mai, p.y_E_CH_NO,
                               p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat]) = catabolic_fluxes(ureserve, Aturnover, r)
    (J1[:C,:rej], J1[:N,:rej], J1[:E,:cat]) =
    synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)
    return (0.0, r)
end
function catabolism!(s::S, u::AbstractStateCNE, t::Float64)::NTuple{2,Float64} where S
    p = s.params
    J1 = s.J1
    Aturnover = (p.k_EC, p.k_EN, p.k_E) .* s.A
    ureserve = (u.C, u.N, u.E)
    m = ureserve ./ u.V
    r = find_rate(t, s.rates, (m, Aturnover, p.j_E_mai, p.y_E_CH_NO,
                               p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat], J1[:EE,:cat]) = catabolic_fluxes(ureserve, Aturnover, r)
    (J1[:C,:rej], J1[:N,:rej], J1[:CN,:cat]) =
    synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)

    J1[:E,:cat] = J1[:EE,:cat] + J1[:CN,:cat] # Total catabolic flux
    θE = J1[:EE,:cat]/J1[:E,:cat] # Proportion of general reserve flux in total catabolic flux
    return (θE, r)
end


"""
Dissipation for any reserve.
Growth, maturity and maintenence are grouped as dissipative processes.
"""
function dissipation!(s::S, u::AbstractState, θE::Float64, r::Float64)::Void where S
    growth!(s, u, θE, r)
    J_E_rep_mai = maturity!(s, u, θE)
    maintenence!(s, u, θE, J_E_rep_mai)
    return nothing
end

function growth!(s::S, u::AbstractState, θE::Float64, r::Float64)::Void where S
    p = s.params; J = s.J; J1 = s.J1;
    y_E_V = 1/p.y_V_E
    y_loss = 1 - p.y_P_V - p.y_V_E

    J[:V,:gro] = r * u.V # Growth flux
    J[:P,:gro] = J[:V,:gro] * p.y_P_V
    drain = -y_E_V * r * u.V
    reserve_drain!(J, u, :gro, drain, θE, p)
    reserve_loss!(J1, J, u, :gro, 1.0)
    J1[:E,:los] -= J[:V,:gro] + J[:P,:gro]
    return nothing
end

# """
#     maturity!(s, u, θE::Float64)
# Calculates reserve drain due to maturity maintenance.
# Stores in M state variable if it exists.
# """
@traitfn function maturity!{S,X; !StateHasM{X}}(s::S, u::X, θE::Float64)
    p = s.params; J = s.J; J1 = s.J1;
    # Only run for structures that take part in reproduction.
    if p.κrep > 0.0
        # TODO: why does rep maintenance stop increasing at M_Vrep?
        # Is this a half finished reproduction model?
        # TODO: -θE * E_rep_mai here, but also θE * E_rep_mai included
        # below in maintenance. Is this intended?

        # Maturity maintenance costs
        J_E_rep_mai = -p.j_E_rep_mai * min(u.V, p.M_Vrep)
        drain = -(p.κrep * J1[:E,:cat] + J_E_rep_mai)
        reserve_drain!(J, u, :rep, drain, θE, p)
        reserve_loss!(J1, J, u, :rep, 1.0)
    else
        J_E_rep_mai = 0.0
    end
    return J_E_rep_mai
end
@traitfn function maturity!{S,X; StateHasM{X}}(s::S, u::X, θE::Float64)
    p = s.params; J = s.J; J1 = s.J1;
    # Only run for structures that take part in reproduction.
    if p.κrep > 0.0
        # Maturity maintenance costs
        J_E_rep_mai = -p.j_E_rep_mai * min(u.V, p.M_Vrep)
        J[:M,:gro] = p.κrep * J1[:E,:cat]
        drain = -(J[:M,:gro] + J_E_rep_mai)
        reserve_drain!(J, u, :rep, drain, θE, p)
        reserve_loss!(J1, J, u, :rep, 1.0)
    else
        J_E_rep_mai = 0.0
    end
    return J_E_rep_mai
end

"""
    maturity!(s, u, θE::Float64)
Calculates reserve drain due to maintenance.
Secretes product state P
"""
function maintenence!(s::S, u::AbstractState, θE::Float64, J_E_rep_mai::Float64)::Void where S
    p = s.params; J = s.J; J1 = s.J1;
    drain = -p.j_E_mai * u.V + J_E_rep_mai #(J_E_rep_mai here cancels out with reproduction)
    reserve_drain!(J, u, :mai, drain, θE, p)
    reserve_loss!(J1, J, u, :mai, 1.0)

    J[:P,:mai] = u.V * p.j_P_mai
    J1[:E,:los] -= J[:P,:mai]

    return nothing
end


"""
    reuse_rejected_reserve!(s::S1, sn::S2) where {S1,S2}
Reallocate state rejected from synthesizing units.
TODO add a 1 structure method
"""
function reuse_rejected_reserve!(s::S1, sn::S2) where {S1,S2}
    p = s.params; J = s.J; J1 = s.J1;
    # This decision seems arbitrary.
    # How would growth dynamics change if both rejected C and N were translocated?
    if s.name == :shoot
        J[:N,:rej] = -p.κEN * J1[:N,:rej]
        J[:C,:rej] = (p.κEC - 1) * J1[:C,:rej]
    else
        J[:C,:rej] = -p.κEC * J1[:C,:rej]
        J[:N,:rej] = (p.κEN - 1) * J1[:N,:rej]
    end

    return nothing
end


"""
    translocate!(s, sn, u::AbstractState)

Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent structures. 
This function is identical both directions, so sn represents
whichever is not the current structure. Will not run with less than 2 structures.

FIXME this will be broken for structure > 2
"""
function translocate!(s::S1, sn::S2, u::AbstractState)::Void where {S1,S2}
    return nothing
end
function translocate!(s::S1, sn::S2, u::AbstractStateE)::Void where {S1,S2}
    κtra = calc_κtra(s.params)
    κtraT = calc_κtra(sn.params)
    # Translocation of C, N and general reserves between structures
    Cgain = κtraT * sn.J1[:E,:cat]
    drain = -κtra * s.J1[:E,:cat]
    s.J[:E,:tra] += Cdrain + Cgain
    s.J1[:E,:los] -= gain
    return nothing
end
function translocate!(s::S1, sn::S2, u::AbstractStateCN)::Void where {S1,S2}
    κtra = calc_κtra(s.params)
    κtraT = calc_κtra(sn.params)

    # Translocation of C, N and general reserves between structures
    tra_C = -κtra * s.J1[:C,:cat]
    tra_N = -κtra * s.J1[:N,:cat]

    s.J[:C,:tra] = tra_C * 1/s.params.y_E_CH_NO
    Cdrain = -κtra * s.J1[:C,:cat]
    Cgain = s.params.y_E_ET * κtraT * sn.J1[:C,:cat]
    s.J[:C,:tra] = Cdrain + Cgain

    s.J[:N,:tra] = tra_N * 1/s.params.y_E_EN
    Ndrain = -κtra * s.J1[:N,:cat]
    Ngain = s.params.y_E_ET * κtraT * sn.J1[:N,:cat]
    s.J[:N,:tra] = Ndrain + Ngain
    return nothing
end
function translocate!(s::S1, sn::S2, u::AbstractStateCNE)::Void where {S1,S2}
    κtra = calc_κtra(s.params)
    κtraT = calc_κtra(sn.params)

    θE = s.J1[:EE,:cat]/s.J1[:E,:cat]

    drain = -κtra * s.J1[:E,:cat]
    reserve_drain!(s.J, u, :tra, drain, θE, s.params)
    reserve_loss!(s.J1, s.J, u, :tra, 1.0)

    gain = s.params.y_E_ET * κtraT * sn.J1[:E,:cat]
    # Final general reserve status
    s.J[:E,:tra] += gain
    s.J1[:E,:los] -= gain
    return nothing
end



"""
    reserve_drain!(J, u, col, drain, θE, params)
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
function reserve_drain!(J::Flux, u::AbstractState, col::Symbol, drain::Float64,
                        θE::Float64, params::P)::Void where P
    return nothing
end
function reserve_drain!(J::Flux, u::AbstractStateE, col::Symbol, drain::Float64,
                        θE::Float64, params::P)::Void where P
    J[:E,col] = drain * θE
    return nothing
end
function reserve_drain!(J::Flux, u::AbstractStateCN, col::Symbol, drain::Float64,
                        θE::Float64, params::P)::Void where P
    J[:C,col] = drain/params.y_E_CH_NO
    J[:N,col] = drain/params.y_E_EN
    return nothing
end
function reserve_drain!(J::Flux, u::AbstractStateCNE, col::Symbol, drain::Float64,
                        θE::Float64, params::P)::Void where P
    J_CN = drain * (1.0 - θE)
    J[:C,col] = J_CN/params.y_E_CH_NO
    J[:N,col] = J_CN/params.y_E_EN
    J[:E,col] = drain * θE
    return nothing
end


"""
    reserve_loss!(J1, J, u, col, drain, θE, params)
Generalised reserve loss to track carbon. 
"""
function reserve_loss!(J1::Flux, J::Flux, u::AbstractStateE, col::Symbol, θloss)::Void
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end
function reserve_loss!(J1::Flux, J::Flux, u::AbstractStateCN, col::Symbol, θloss)::Void
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    return nothing
end
function reserve_loss!(J1::Flux, J::Flux, u::AbstractStateCNE, col::Symbol, θloss)::Void
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end


"""
    calc_κtra(params::P) where P
κtra is the difference between κsoma and κrep
"""
function calc_κtra(params::P) where P
    1.0 - params.κsoma - params.κrep
end


# J: Flux matrix diagram.
# Rows: state.
# Columns: transformations
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ assS       │ groS       │ maiS       │ repS       │ rejS       │ traS       ┃ assR       │ groR       │ maiR       │ repR │ rejR       │ traR       ┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃    ┃ JSS SubArray                                                                ┃ JRR SubArray                                                          ┃
# ┃    ┃                                                                             ┃                                                                       ┃
# ┃PS  ┃ 0          │ J_PS,groS  │ J_PS,maiS  │ 0          │ 0          │ 0          ┃ 0          │ J_PR,groR  │ J_PR,maiR  │ 0    │ 0          │ 0          ┃
# ┃VS  ┃ 0          │ J_VS,groS  │ 0          │ 0          │ 0          │ 0          ┃ 0          │ J_VR,groR  │ 0          │ 0    │ 0          │ 0          ┃
# ┃RS  ┃ 0          │ 0          │ 0          │ J_MS,groS  │ 0          │ 0          ┃ 0          │ 0          │ 0          │ 0    │ 0          │ 0          ┃
# ┃ECS ┃ J_ECS,assS │ J_ECS,groS │ J_ECS,maiS │ J_ECS,repS │ J_ECS,rejS │ J_ECS,traS ┃ J_ECR,assR │ J_ECR,groR │ J_ECR,maiR │ 0    │ J_ECS,rejR │ J_ECR,traR ┃
# ┃ENS ┃ J_ENS,assS │ J_ENS,groS │ J_ENS,maiS │ J_ENS,repS │ J_ENS,rejS │ J_ENS,traS ┃ J_ENR,assR │ J_ENR,groR │ J_ENR,maiR │ 0    │ J_ENS,rejR │ J_ENR,traR ┃
# ┃ES  ┃ J_ES,assS  │ J_ES,groS  │ J_ES,maiS  │ J_ES,repS  │ 0          │ J_ES,traS  ┃ J_ER,assR  │ J_ER,groR  │ J_ER,maiR  │ 0    │ 0          │ J_ER,traR  ┃
# ┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# J1: Catabolic flux matrix diagrams.
# ┏━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃     ┃ catS       │ rejS      │ losS      ┃
# ┣━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃EES  ┃ J_EES,catS │ 0         │ 0         ┃
# ┃CNS  ┃ J_CNS,catS │ 0         │ 0         ┃
# ┃CS   ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃
# ┃NS   ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃
# ┃ES   ┃ J_ES,catS  │ 0         │ 0         ┃
# ┃EES  ┃ J_EER,catR │ 0         │ 0         ┃
# ┃NS   ┃ J_CNR,catR │ 0         │ 0         ┃
# ┃CS   ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
# ┃ENS  ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
# ┃ES   ┃ J_ER,catR  │ 0         │ 0         ┃
# ┗━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# M: State vector diagram.
# ┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                      ┃ MS SubArray            ┃ MR SubArray            ┃
# ┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
