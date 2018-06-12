# import BiophysicalModels.runmodel!

"""
    runmodel!(du, settings, u, t)
DEBSettings method for BiophysicalModels.jl api.
Applies environment and runs the DEB model.
"""
function runmodel!(du, u, p, t::Number)
    organism = p.nodes[1]
    apply(setvars!, organism.nodes, t)
    offset_apply!(setstate!, organism.nodes, u, 0)
    apply(apply_environment!, organism.nodes, p.environment, t)
    apply(zero_flux!, organism.nodes)

    debmodel!(organism, t)
    offset_apply!(sum_flux!, du, organism.nodes, 0)

    return nothing
end

function diff(du, u, p, t::Number)
    organism = p.nodes[1]
    debmodel!(organism, t)
    offset_apply!(sum_flux!, du, organism.nodes)
    du
end

setvars!(o, t) = o.vars = o.varsrecord[t] 

zero_flux!(o) = begin
    fill!(o.J, zero(J[1,1]))
    fill!(o.J1, zero(J1[1,1]))
end

setstate!(o, u, offset::Int) = begin
    for i in 1:length(o.state) 
        o.state[i] = u[i+offset]
    end
    offset + length(o.state)
end

sum_flux!(du, o, offset::Int) = begin
    for i in 1:size(o.J, 1) 
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
    apply(translocation!, organs, swapped)
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
Catabolism for E, C and N, or C, N and E reserves.

Does not finalise flux in J - operates only on J1 (intermediate storage)
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
    maturity!(o.params.maturity, o, u)
    maintenence!(o, u)
    product!(o, u)
    return nothing
end

"""
    growth!(o, u::AbstractState, θE, r)
Allocates reserves to growth.
"""
function growth!(o, u)
    v = o.vars; p = o.params; J = o.J; J1 = o.J1;
    grow = v.rate * u.V
    J[:V,:gro] = grow 
    drain = -(1/p.y_V_E) * grow 
    loss = (1/p.y_V_E - 1) * v.rate * u.V
    reserve_drain!(o, :gro, drain, v.θE)
    reserve_loss!(o, loss)
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
    drain = -(f.κrep * J1[:E,:cat] + -f.j_E_rep_mai * min(u.V, f.M_Vrep))
    reserve_drain!(o, :rep, drain, v.θE)
    reserve_loss!(o, -drain)
    return nothing
end

@traitfn function maturity!{X; StateHasM{X}}(f::Maturity, o, u::X)
    v = o.vars; p = o.params; J = o.J; J1 = o.J1;
    J[:M,:gro] = f.κrep * J1[:E,:cat]
    maint = -f.j_E_rep_mai * u.V
    drain = -J[:M,:gro] + maint # min(u.V, f.M_Vrep))
    reserve_drain!(o, :rep, drain, v.θE)
    reserve_loss!(o, -maint)
    return nothing
end

"""
    maintenance!(o, u, θE)
Allocates reserve drain due to maintenance.
"""
function maintenence!(o, u)
    drain = -o.params.j_E_mai * u.V
    reserve_drain!(o, :mai, drain, o.vars.θE)
    reserve_loss!(o, -drain) # all maintenance is loss
    nothing
end

"""
    product!(o, u::AbstractState, θE, r)
Allocates waste products from growth and maintenance.
"""
function product!(o, u)
    o.J[:P,:gro] = o.J[:V,:gro] * o.params.y_P_V
    o.J[:P,:mai] = u.V * o.params.j_P_mai
    # undo the reserve loss from growth: it went to product
    loss_correction = -(o.J[:P,:gro] + o.J[:P,:mai])
    reserve_loss!(o, loss_correction)
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
function translocation!(o, on)
    reuse_rejected!(o, on)
    translocate!(o, on)
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
function translocate!(o, on)
    p = o.params; J = o.J; J1 = o.J1;
    trans = -κtra(p) * J1[:E,:cat] # fraction of catabolised reserves translocated
    reserve_drain!(o, :tra, trans, o.vars.θE) # drain from translocation
    on.J[:E,:tra] += p.y_E_ET * trans # translocation to other organ
    # rejected reserves are translocated and used in assimilation.
    return nothing
end

"""
    reuse_rejected!(o, on)
Reallocate state rejected from synthesizing units.

TODO add a 1-organs method
Also how does this interact with assimilation?
"""
function reuse_rejected!(o, on)
    p = o.params
    # rejected reserves are translocated and used in assimilation.
    o.J[:C,:rej] = -o.J1[:C,:rej]
    o.J[:N,:rej] = -o.J1[:N,:rej]
    on.J[:C,:tra] = p.y_EC_ECT * o.J1[:C,:rej]
    on.J[:N,:tra] = p.y_EN_ENT * o.J1[:N,:rej]
    o.J1[:C,:los] += (1 - p.y_EC_ECT) * o.J1[:C,:rej]
    o.J1[:N,:los] += (1 - p.y_EN_ENT) * o.J1[:N,:rej]
    return nothing
end

"""
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
function reserve_drain!(o, col, drain, θE)
    J_CN = drain * (1.0 - θE) # fraction on drain from C and N reserves
    o.J[:C,col] = J_CN/o.params.y_E_CH_NO
    o.J[:N,col] = J_CN/o.params.y_E_EN
    o.J[:E,col] = drain * θE
    return nothing
end

"""
Generalised reserve loss to track carbon. 
"""
function reserve_loss!(o, loss)
    o.J1[:C,:los] += loss/o.params.y_E_CH_NO
    o.J1[:N,:los] += loss/o.params.y_E_EN
    return nothing
end

"""
    κtra(params::P) where P
κtra is the difference paramsbetween κsoma and κrep
"""
κtra(params) = 1.0 - params.κsoma - κrep(params)
κrep(params) = :κrep in fieldnames(params.maturity) ? params.maturity.κrep : 0.0

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
