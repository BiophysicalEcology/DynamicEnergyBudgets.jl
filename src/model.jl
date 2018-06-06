# import BiophysicalModels.runmodel!

"""
    runmodel!(du, settings, u, t)
DEBSettings method for BiophysicalModels.jl api.
Applies environment and runs the DEB model.
"""
function runmodel!(du, u, p, t::Number)
    organism = p.nodes[1]
    apply(setvars!, organism.nodes, t)
    offset_apply!(setstate!, u, organism.nodes)
    apply(apply_environment!, organism.nodes, p.environment, t)

    debmodel!(organism, t)
    offset_apply!(sum_flux!, du, organism.nodes)

    return nothing
end

function diff(du, u, p, t::Number)
    organism = p.nodes[1]
    debmodel!(organism, t)
    offset_apply!(sum_flux!, du, organism.nodes)
    du
end

setvars!(o, t) = o.vars = o.varsrecord[t] 

setstate!(u, offset::Int, o) = begin
    for i in 1:length(o.state) 
        o.state[i] = u[i+offset]
    end
    offset + length(o.state)
end

sum_flux!(du, offset::Int, o) = begin
    for i in 1:length(o.state) 
        du[i+offset] = sum(o.J[i,:])
    end
    offset + length(o.state)
end


"""
    deb_model!(settings, t)
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
function debmodel!(organism, t::Number)
    organs = organism.nodes
    swapped = (Base.tail(organs)..., organs[1])

    apply(metabolism!, organs, t)
    # FIXME what about rejected reserves with 1 organ?
    length(organs) > 1 && apply(translocation!, organs, swapped)
    apply(assimilation!, organs, swapped)
    return nothing
end

"""
    metabolism!(s, t::Number)
Metabolism is an identical process for all organs, with potentially
different parameters or area and rate functions.
"""
function metabolism!(o, t::Number)
    o.vars.scale = scaling(o.params.scaling, o.state.V)
    catabolism!(o, o.state, t)
    dissipation!(o, o.state)
    feedback!(o, o.params.feedback, o.state)
    return nothing
end

"""
    translocation!(s, on)
Some rejected reserve is translocated.
But this grouping of functions is somewhat unsatisfying.
"""
function translocation!(o, on)
    reuse_rejected_reserve!(o, on)
    translocate!(o, on, o.state)
    return nothing
end

"""
Catabolism for E, C and N, or C, N and E reserves.
"""
function catabolism!(o, u::AbstractStateCNE, t::Number)
    p = o.params; v = o.vars; J1 = o.J1
    scaledturnover = (p.k_EC, p.k_EN, p.k_E) .* v.scale
    ureserve = (u.C, u.N, u.E)
    m = ureserve ./ u.V
    v.rate = find_rate(v, (m, scaledturnover, p.j_E_mai, p.y_E_CH_NO, p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat], J1[:EE,:cat]) = 
        catabolic_fluxes(ureserve, scaledturnover, v.rate)
    (J1[:C,:rej], J1[:N,:rej], J1[:CN,:cat]) =
        synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)
    J1[:E,:cat] = J1[:EE,:cat] + J1[:CN,:cat] # Total catabolic flux
    v.θE = J1[:EE,:cat]/J1[:E,:cat] # Proportion of general reserve flux in total catabolic flux
    return nothing
end

"""
    dissipation!(s, u::AbstractState, θE, r)
Dissipation for any reserve.
Growth, maturity and maintenence are grouped as dissipative processes.
"""
function dissipation!(o, u)
    growth!(o, u)
    maturity!(o.maturity, o, u)
    maintenence!(o, u)
    production!(o, u)
    return nothing
end

"""
    growth!(o, u::AbstractState, θE, r)
Allocates reserves to growth.
"""
function growth!(o, u)
    v = o.vars; p = o.params; J = o.J; J1 = o.J1;
    y_E_V = 1/p.y_V_E
    J[:V,:gro] = v.rate * u.V # Growth flux
    drain = -y_E_V * v.rate * u.V
    reserve_drain!(J, u, :gro, drain, v.θE, p)
    reserve_loss!(J1, J, u, :gro, 1.0)
    J1[:E,:los] -= J[:V,:gro]
    return nothing
end

"""
    maturity!(o, u, θE)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
function maturity!(f, o, u) end

@traitfn function maturity!{X; !StateHasM{X}}(f::Maturity, o, u::X)
    v = o.vars; p = o.params; J = o.J; J1 = o.J1;
    # TODO: why does rep maintenance stop increasing at M_Vrep?
    # Is this a half finished reproduction model?
    # TODO: -θE * E_rep_mai here, but also θE * E_rep_mai included
    # below in maintenance. Is this intended?
    drain = -(p.κrep * J1[:E,:cat] + -p.j_E_rep_mai * min(u.V, p.M_Vrep))
    reserve_drain!(J, u, :rep, drain, v.θE, p)
    reserve_loss!(J1, J, u, :rep, 1.0)
    return nothing
end
@traitfn function maturity!{X; StateHasM{X}}(f::Maturity, u::X)
    v = o.vars; p = o.params; J = o.J; J1 = o.J1;
    J[:M,:gro] = f.κrep * J1[:E,:cat]
    drain = -(J[:M,:gro] + -f.j_E_rep_mai * min(u.V, f.M_Vrep))
    reserve_drain!(J, u, :rep, drain, v.θE, p)
    reserve_loss!(J1, J, u, :rep, 1.0)
    return nothing
end

"""
    maintenance!(o, u, θE)
Allocates reserve drain due to maintenance.
"""
function maintenence!(o, u)
    v = o.vars; p = o.params 
    drain = -p.j_E_mai * u.V
    reserve_drain!(o.J, u, :mai, drain, v.θE, p)
    reserve_loss!(o.J1, o.J, u, :mai, 1.0)
    return nothing
end

"""
    production!(o, u::AbstractState, θE, r)
Allocates waste products from growth and maintenance.
"""
function production!(o, u)
    p = o.params; J = o.J; J1 = o.J1;
    J[:P,:gro] = J[:V,:gro] * p.y_P_V
    J[:P,:mai] = u.V * p.j_P_mai
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
    p = o.params; J = o.J; J1 = o.J1;
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
translocate!(o, on, u) = nothing

function translocate!(o, on, u::AbstractStateCNE)
    p = o.params; J = o.J; J1 = o.J1;

    drain = -calc_κtra(p) * J1[:E,:cat]
    reserve_drain!(J, u, :tra, drain, o.vars.θE, p)
    reserve_loss!(J1, J, u, :tra, 1.0)

    gain = p.y_E_ET * calc_κtra(on.params) * on.J1[:E,:cat]
    J[:E,:tra] += gain
    J1[:E,:los] -= gain
    return nothing
end

"""
    reserve_drain!(J, u, col, drain, θE, params)
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
function reserve_drain!(J, u::AbstractStateCNE, col, drain, θE, params)
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
function reserve_loss!(J1, J, u::AbstractStateCNE, col, θloss)
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end

"""
    calc_κtra(params::P) where P
κtra is the difference between κsoma and κrep
"""
calc_κtra(params) = 1.0 - params.κsoma - κrep(params)

κrep(params) = typeof(params.maturity) == nothing ? 0.0 : params.maturity.κrep

"""
    find_rate(t, args)
Calculate rate formula. TODO: use Roots.jl for this
"""
function find_rate(v, args::Tuple{NTuple{N},NTuple{N},Vararg}) where {N}
    local f = rate_formula
    local found = false

    # Find the largest possible rate window.
    bounds = (-1.0oneunit(v.rate), 2.0oneunit(v.rate)) #rate_window(args...) 
    find_zero(x -> rate_formula(x, args...), bounds, Secant(); xtol=BI_XTOL)
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
