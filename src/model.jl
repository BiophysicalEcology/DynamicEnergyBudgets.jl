import BiophysicalModels.runmodel!

"""
    runmodel!(du, settings, u, t)
DEBSettings method for BiophysicalModels.jl api.
Applies environment and runs the DEB model.
"""
function runmodel!(du, u, scenario::Scenario, t::Number)
    organism = scenario.nodes[1]
    apply(apply_environment!, organism, scenario.environment, t)

    deb_model!(organism, t)
    sum_flux!(du, organism) # FIXME du is all organs

    return nothing
end

"""
    deb_model!(settings, t)
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
function deb_model!(organism, t::Number)
    organs = organism.nodes
    swapped = (Base.tail(organs)..., organs[1])

    apply(metabolism!, organs, t)
    # FIXME what about rejected reserves with 1 organ?
    length(organs) > 1 && apply(translocation!, organs, swapped)
    apply(assimilation!, organs, swapped, settings)

    return nothing
end

"""
    metabolism!(s, t::Number)
Metabolism is an identical process for all organs, with potentially
different parameters or area and rate functions.
"""
function metabolism!(o, t::Number)
    s.scaling = scaling(o.scaling, o.u.V)
    (θE, r) = catabolism!(o, o.values, t)
    dissipation!(o, o.values, θE, r)
    feedback!(o, o.feedback, o.values, r)
    return nothing
end

"""
    translocation!(s, on)
Some rejected reserve is translocated.
But this grouping of functions is somewhat unsatisfying.
"""
function translocation!(o, on)
    reuse_rejected_reserve!(o, on)
    translocate!(o, on, o.values)
    return nothing
end

"""
Catabolism for E, C and N, or C, N and E reserves.
"""
function catabolism!(o, u::AbstractStateE, t::Number)
    p = o.params
    J1 = o.values.J1
    Aturnover = p.k * o.A
    m = u.E / u.V
    r = find_rate(t, o.rates, (m, Aturnover, p.j_E_mai, p.y_V_E, p.κsoma))
    J1[:E,:cat] = catabolic_fluxes(ureserve, Aturnover, r)
    return (0.0, r)
end
function catabolism!(o, u::AbstractStateCN, t::Number)
    p = o.params
    J1 = o.values.J1
    Aturnover = p.k .* o.A
    ureserve = (u.C, u.N)
    m = ureserve ./ u.V
    r = find_rate(t, o.rates, (m, Aturnover, p.j_E_mai, p.y_E_CH_NO,
                               p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat]) = catabolic_fluxes(ureserve, Aturnover, r)
    (J1[:C,:rej], J1[:N,:rej], J1[:E,:cat]) =
        synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)
    return (0.0, r)
end
function catabolism!(o, u::AbstractStateCNE, t::Number)
    p = o.params
    J1 = o.values.J1
    Aturnover = p.k .* o.A
    ureserve = (u.C, u.N, u.E)
    m = ureserve ./ u.V
    r = find_rate(t, o.rates, (m, Aturnover, p.j_E_mai, p.y_E_CH_NO,
                               p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat], J1[:EE,:cat]) = catabolic_fluxes(ureserve, Aturnover, r)
    (J1[:C,:rej], J1[:N,:rej], J1[:CN,:cat]) =
        synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)

    J1[:E,:cat] = J1[:EE,:cat] + J1[:CN,:cat] # Total catabolic flux
    θE = J1[:EE,:cat]/J1[:E,:cat] # Proportion of general reserve flux in total catabolic flux
    return (θE, r)
end

"""
    dissipation!(s, u::AbstractState, θE, r)
Dissipation for any reserve.
Growth, maturity and maintenence are grouped as dissipative processes.
"""
function dissipation!(o, u::AbstractState, θE, r)
    growth!(o, u, θE, r)
    maturity!(o, u, θE)
    maintenence!(o, u, θE)
    production!(o, u, θE)
    return nothing
end

"""
    growth!(o, u::AbstractState, θE, r)

Allocates reserves to growth.
"""
function growth!(o, u::AbstractState, θE, r)
    p = o.params; J = o.values.J; J1 = o.values.J1;
    y_E_V = 1/p.y_V_E

    J[:V,:gro] = r * u.V # Growth flux
    drain = -y_E_V * r * u.V
    reserve_drain!(J, u, :gro, drain, θE, p)
    reserve_loss!(J1, J, u, :gro, 1.0)
    J1[:E,:los] -= J[:V,:gro]
    return nothing
end

"""
    maturity!(o, u, θE)

Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
function maturity!(o, u, θE) end
@traitfn function maturity!{X; !StateHasM{X}}(o, u::X, θE)
    p = o.params; J = o.values.J; J1 = o.values.J1;
    # Only run for organs that take part in reproduction.
    if p.κrep > 0.0
        # TODO: why does rep maintenance stop increasing at M_Vrep?
        # Is this a half finished reproduction model?
        # TODO: -θE * E_rep_mai here, but also θE * E_rep_mai included
        # below in maintenance. Is this intended?
        # Maturity maintenance costs
        drain = -(p.κrep * J1[:E,:cat] + -p.j_E_rep_mai * min(u.V, p.M_Vrep))
        reserve_drain!(J, u, :rep, drain, θE, p)
        reserve_loss!(J1, J, u, :rep, 1.0)
    end
    return nothing
end
@traitfn function maturity!{X; StateHasM{X}}(o, u::X, θE)
    p = o.params; J = o.values.J; J1 = o.values.J1;
    # Only run for organs that take part in reproduction.
    # TODO: no logic with parameter values.
    if p.κrep > 0.0
        # Maturity maintenance costs
        J[:M,:gro] = p.κrep * J1[:E,:cat]
        drain = -(J[:M,:gro] + -p.j_E_rep_mai * min(u.V, p.M_Vrep))
        reserve_drain!(J, u, :rep, drain, θE, p)
        reserve_loss!(J1, J, u, :rep, 1.0)
    end
    return nothing
end

"""
    maintenance!(o, u, θE)
Allocates reserve drain due to maintenance.
"""
function maintenence!(o, u::AbstractState, θE)
    p = o.params 
    drain = -p.j_E_mai * u.V # Not temp related. Why?
    reserve_drain!(o.J, u, :mai, drain, θE, p)
    reserve_loss!(o.J1, o.J, u, :mai, 1.0)
    return nothing
end

"""
    production!(o, u::AbstractState, θE, r)

Allocates waste products from growth and maintenance.
"""
function production!(o, u::AbstractState, θE, r)
    p = o.params; J = o.values.J; J1 = o.values.J1;
    J[:P,:gro] = J[:V,:gro] * p.y_P_V
    J[:P,:mai] = u.V * p.j_P_mai # Maintenance product is not temp related. Why?
    J1[:E,:los] -= J[:P,:gro] + J[:P,:mai]
    return nothing
end

"""
    reuse_rejected_reserve!(o, on)
Reallocate state rejected from synthesizing units.

TODO add a 1-organs method
Also how does this interact with assimilation?
"""
function reuse_rejected_reserve!(o, on)
    p = o.params; J = o.values.J; J1 = o.values.J1;
    # This decision seems arbitrary.
    # How would growth dynamics change if both rejected C and N were translocated?
    # FIXME: balance this thing
    if typeof(p.assimilation) <: CarbonAssimilation
        J[:N,:rej] = -p.κEN * J1[:N,:rej]
        J[:C,:rej] = (p.κEC - 1) * J1[:C,:rej]
    elseif typeof(p.assimilation) <: NitrogenAssimilation
        J[:C,:rej] = -p.κEC * J1[:C,:rej]
        J[:N,:rej] = (p.κEN - 1) * J1[:N,:rej]
    end
    return nothing
end

"""
    translocate!(o, on, u::AbstractState)

Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs. 
This function is identical both directiono, so on represents
whichever is not the current organs. Will not run with less than 2 organs.

FIXME this will be broken for organs > 2
"""
function translocate!(o, on, u::AbstractState)
    return nothing
end
function translocate!(o, on, u::AbstractStateE)
    κtra = calc_κtra(o.params)
    κtraT = calc_κtra(on.params)
    # Translocation of C, N and general reserves between organs
    Cgain = κtraT * on.values.J1[:E,:cat]
    drain = -κtra * o.values.J1[:E,:cat]
    o.values.J[:E,:tra] += Cdrain + Cgain
    o.values.J1[:E,:los] -= gain
    return nothing
end
function translocate!(o, on, u::AbstractStateCN)
    J = o.values.J; J1 = o.values.J1;
    Jn = on.values.J; J1n = on.values.J1;
    κtra = calc_κtra(o.params)
    κtraT = calc_κtra(on.params)

    # Translocation of C, N and general reserves between organs
    tra_C = -κtra * J1[:C,:cat]
    tra_N = -κtra * J1[:N,:cat]

    J[:C,:tra] = tra_C * 1/o.paramo.y_E_CH_NO
    Cdrain = -κtra * J1[:C,:cat]
    Cgain = o.params.y_E_ET * κtraT * J1n[:C,:cat]
    J[:C,:tra] = Cdrain + Cgain

    J[:N,:tra] = tra_N * 1/o.params.y_E_EN
    Ndrain = -κtra * J1[:N,:cat]
    Ngain = o.params.y_E_ET * κtraT * J1n[:N,:cat]
    J[:N,:tra] = Ndrain + Ngain
    return nothing
end
function translocate!(o, on, u::AbstractStateCNE)
    p = o.params, J = o.values.J; J1 = o.values.J1;

    θE = J1[:EE,:cat]/J1[:E,:cat]

    drain = -calc_κtra(p) * J1[:E,:cat]
    reserve_drain!(J, u, :tra, drain, θE, p)
    reserve_loss!(J1, J, u, :tra, 1.0)

    gain = p.y_E_ET * calc_κtra(on.params) * on.values.J1[:E,:cat]
    J[:E,:tra] += gain
    J1[:E,:los] -= gain
    return nothing
end

"""
    reserve_drain!(J, u, col, drain, θE, params)
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
function reserve_drain!(J, u::AbstractState, col::Symbol, drain, θE, params)
    return nothing
end
function reserve_drain!(J, u::AbstractStateE, col::Symbol, drain, θE, params)
    J[:E,col] = drain * θE
    return nothing
end
function reserve_drain!(J, u::AbstractStateCN, col::Symbol, drain, θE, params)
    J[:C,col] = drain/params.y_E_CH_NO
    J[:N,col] = drain/params.y_E_EN
    return nothing
end
function reserve_drain!(J, u::AbstractStateCNE, col::Symbol, drain, θE, params)
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
function reserve_loss!(J1, J, u::AbstractStateE, col::Symbol, θloss)
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end
function reserve_loss!(J1, J, u::AbstractStateCN, col::Symbol, θloss)
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    return nothing
end
function reserve_loss!(J1, J, u::AbstractStateCNE, col::Symbol, θloss)
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end

"""
    calc_κtra(params::P) where P
κtra is the difference between κsoma and κrep
"""
calc_κtra(params) = 1.0 - params.κsoma - params.κrep

"""
    find_rate(t, rates::Array{Float64}, args)
Calculate rate formula.

TODO: use Roots.jl for this
"""
function find_rate(t, rates, args::Tuple{NTuple{N},NTuple{N},Vararg}) where {N}
    local f = rate_formula
    local found = false

    # Find the largest possible rate window.
    x0, x1 = (-1.0u"mol/mol*d^-1", 2.0u"mol/mol*d^-1") #rate_window(args...) 

    # bounds = (-1.0u"mol/mol*d^-1", 2.0u"mol/mol*d^-1") #rate_window(args...) 
    find_zero(x -> rate_formula(x, args...), bounds, FalsePosition(); xtol=BI_XTOL)

    # Save rate for plotting.
    save_rate!(t, rates, x)
    return x
end

function save_rate!(t, rates, newrate)
    t_floor = floor(Int64, t)
    if t_floor == t
        rates[t_floor + 1] = ustrip(newrate)
    end
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
