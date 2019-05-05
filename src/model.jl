"""
Run a DEB organism model from DiffEq.
"""
(organism::Plant)(du, u, t::Number, organs) = begin
    dead(organism) && return

    # Make sure the parameters don't any physical laws
    check_params(organs)

    # Split state into separate organs
    ux = split_state(organs, u)

    # Set up variables for this timestep and the current state
    apply(update_tstep!, organs, t)
    apply(zero_flux!, organs)
    apply(update_height!, organs, ux)
    apply(update_shape!, organs, ux)
    apply_environment!(organism, organs, ux, t)

    # Run the model, tag the organism as dead if it breaks.
    if !debmodel!(organs, ux, environment(organism)) 
        du .= zero(eltype(du)) 
        set_dead!(organism, true)
        return
    end

    # Sum the flux matrix to the state change vector
    sum_flux!(du, organs)
    nothing
end

"""
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation methods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
debmodel!(organs::Tuple, u::Tuple, env) = begin
    # Quit if it dies 
    false in metabolism!(organs, u) && return false
    translocation!(organs)
    assimilation!(organs, u)
    true
end

"""
Metabolism is an identical process for all organs, with potentially
different parameters or area and rate functions.
"""
metabolism!(organs::Tuple, u) = apply(metabolism!, organs, u)
metabolism!(o::AbstractOrgan, u) = begin
    # Quit if it dies 
    catabolism!(o, u) || return false
    growth!(o, u)
    maturity!(o, u)
    maintenence!(o, u)
    resorption!(o, u)
    return true
end

"""
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

"""
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
@inline reserve_drain!(o::AbstractOrgan, col, drain) = reserve_drain!(has_reserves(o), o, col, drain)
@inline reserve_drain!(::HasCNE, o, col, drain) = begin
    ΘE, J = θE(o), flux(o)
    J_CN = -drain * (1 - θ) # fraction on drain from C and N reserves
    @inbounds J[:C,col] = J_CN/y_E_EC(o)
    @inbounds J[:N,col] = J_CN/y_E_EN(o)
    @inbounds J[:E,col] = -drain * ΘE 
    nothing
end
@inline reserve_drain!(::HasCN, o, col, drain) = begin
    J = flux(o)
    @inbounds J[:C,col] = -drain/y_E_EC(o)
    @inbounds J[:N,col] = -drain/y_E_EN(o)
    nothing
end

# J: Flux matrix diagram.
# Rows: state.
# Columns: transformations
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ assS       │ groS       │ maiS       │ matS       │ rejS       │ traS       ┃ assR       │ groR       │ maiR       │ ratR │ rejR       │ traR       ┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃    ┃ JSS SubArray                                                                ┃ JRR SubArray                                                          ┃
# ┃    ┃                                                                             ┃                                                                       ┃
# ┃PS  ┃ 0          │ J_PS,groS  │ J_PS,maiS  │ 0          │ 0          │ 0          ┃ 0          │ J_PR,groR  │ J_PR,maiR  │ 0    │ 0          │ 0          ┃
# ┃VS  ┃ 0          │ J_VS,groS  │ 0          │ 0          │ 0          │ 0          ┃ 0          │ J_VR,groR  │ 0          │ 0    │ 0          │ 0          ┃
# ┃RS  ┃ 0          │ 0          │ 0          │ J_MS,groS  │ 0          │ 0          ┃ 0          │ 0          │ 0          │ 0    │ 0          │ 0          ┃
# ┃ECS ┃ J_ECS,assS │ J_ECS,groS │ J_ECS,maiS │ J_ECS,matS │ J_ECS,rejS │ J_ECS,traS ┃ J_ECR,assR │ J_ECR,groR │ J_ECR,maiR │ 0    │ J_ECS,rejR │ J_ECR,traR ┃
# ┃ENS ┃ J_ENS,assS │ J_ENS,groS │ J_ENS,maiS │ J_ENS,matS │ J_ENS,rejS │ J_ENS,traS ┃ J_ENR,assR │ J_ENR,groR │ J_ENR,maiR │ 0    │ J_ENS,rejR │ J_ENR,traR ┃
# ┃ES  ┃ J_ES,assS  │ J_ES,groS  │ J_ES,maiS  │ J_ES,matS  │ 0          │ J_ES,traS  ┃ J_ER,assR  │ J_ER,groR  │ J_ER,maiR  │ 0    │ 0          │ J_ER,traR  ┃
# ┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# J1: Catabolic flux diagram.
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ catS       │ rejS      │ losS mols!┃    ┃ catR       │ rejR      │ losR mols!┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃EES ┃ J_EES,catS │ 0         │ 0         ┃EES ┃ J_EER,catR │ 0         │ 0         ┃
# ┃CNS ┃ J_CNS,catS │ 0         │ 0         ┃NS  ┃ J_CNR,catR │ 0         │ 0         ┃
# ┃CS  ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃CS  ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
# ┃NS  ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃ENS ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
# ┃ES  ┃ J_ES,catS  │ 0         │ 0         ┃ES  ┃ J_ER,catR  │ 0         │ 0         ┃
# ┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# M: State vector diagram.
# ┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
