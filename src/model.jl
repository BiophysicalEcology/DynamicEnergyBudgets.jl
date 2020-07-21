
"""
    metabolism!(organs::Tuple, u::Tuple)
    metabolism!(o::AbstractOrgan, u::AbstractVector)

Metabolism is the same basic process for all organs, 
with potentially different components and parameters.

`catabolism!` determines the current growth rate, flux
from reserve and wether the organism is still alive,

Then growth, maintenence, maturity and resorption update
the flux matrix `o.J` based on catabolised reserve and state `u`.
"""
function metabolism! end
metabolism!(organs::Tuple, u::Tuple) = map(metabolism!, organs, u)
metabolism!(o::AbstractOrgan, u::AbstractVector) = begin
    alive = catabolism!(o, u) 
    alive || return false

    growth!(o, u)
    maintenence!(o, u)
    maturity!(o, u)
    resorption!(o, u)
    return true
end

"""
    debmodel!(organs::Tuple, u::Tuple{AbstractArray}, env)

A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

The method applies metabolism, translocation and assimilation methods 
to all organs. If metabolism fails in any organ, the organism is dead,
and `false` is returned.

`organs` is a tuple of Organ, `u` is a tuple of organ state variable vectors, 
env is the environment component (or `nothing`).
"""
function debmodel!(organs::Tuple, u::Tuple, env)
    # Quit if it dies
    false in metabolism!(organs, u) && return false
    translocation!(organs)
    assimilation!(organs, u)
    return true
end

# These are currently implemented for Plant. In future this will be implemented
# for AbstractOrganism, when julis allows that.
(o::Plant)(du::AbstractVector{<:Unitful.Quantity}, u::AbstractVector{<:Unitful.Quantity}, p, t::Unitful.Quantity) = begin
    DynamicEnergyBudgets.dead(o) && return
    o(du, u, t, update_organs(organs(o), t))
    return
end
# Deal with du, u, p without units
(o::Plant)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Real) = begin
    DynamicEnergyBudgets.dead(o) && return
    # Create unitful state vectors
    du1 = du .* (mol/hr)
    u1 = u .* mol
    t1 = t * hr
    # Reconstruct the model with new parameters
    o1 = reconstruct(o, p, Real)
    # Run the model
    o1(du1, u1, t1, update_organs(organs(o1), t1))
    # Copy the result to the unitless state change vector
    du2 = du1 ./ (mol/hr)
    if eltype(du2) == eltype(du)
        du .= du2
    end
    return
end
# Deal with du and u without units
(o::Plant)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p, t::Real) = begin 
    DynamicEnergyBudgets.dead(o) && return
    # Create unitful state vectors
    du1 = du .* (mol/hr)
    u1 = u .* mol
    t1 = t * hr
    # Run the model
    o(du1, u1, t1, update_organs(organs(o), t1))
    # Copy the result to the unitless state change vector
    du .= du1 ./ (mol/hr)
    return
end
(organism::Plant)(du, u, t::Number, organs::Tuple) = begin
    DynamicEnergyBudgets.dead(organism) && return
    # Make sure the parameters don't break any physical laws
    check_params(organs)

    # Split state into separate organs
    ux = split_state(organs, u)

    # Set up variables for this timestep and the current state
    map(zero_flux!, organs)
    map(update_height!, organs, ux)
    map(update_scaling!, organs, ux)
    try
        apply_environment!(organism, organs, ux, t)
    catch e
        du .= zero(eltype(du))
        set_dead!(organism, true)
        @warn "dead at $t due to error: $e"
        return
    end

    # Run the model, tag the organism as dead if it breaks.
    if !debmodel!(organs, ux, environment(organism))
        du .= zero(eltype(du))
        set_dead!(organism, true)
        @warn "dead at $t due to growth rate"
        return
    end

    # Sum the flux matrix to the state change vector
    sum_flux!(du, organs)
    return
end

#= J: Flux matrix diagram. rows=state, columns=transformations
┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃    ┃ assS       │ groS       │ maiS       │ matS       │ rejS       │ traS       ┃
┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃    ┃ JSS SubArray                                                                ┃
┃    ┃                                                                             ┃
┃PS  ┃ 0          │ J_PS,groS  │ J_PS,maiS  │ 0          │ 0          │ 0          ┃
┃VS  ┃ 0          │ J_VS,groS  │ 0          │ 0          │ 0          │ 0          ┃
┃RS  ┃ 0          │ 0          │ 0          │ J_MS,groS  │ 0          │ 0          ┃
┃ECS ┃ J_ECS,assS │ J_ECS,groS │ J_ECS,maiS │ J_ECS,matS │ J_ECS,rejS │ J_ECS,traS ┃
┃ENS ┃ J_ENS,assS │ J_ENS,groS │ J_ENS,maiS │ J_ENS,matS │ J_ENS,rejS │ J_ENS,traS ┃
┃ES  ┃ J_ES,assS  │ J_ES,groS  │ J_ES,maiS  │ J_ES,matS  │ 0          │ J_ES,traS  ┃
┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃    ┃ assR       │ groR       │ maiR       │ matR       │ rejR       │ traR       ┃
┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃    ┃ JRR SubArray                                                                ┃
┃    ┃                                                                             ┃
┃PS  ┃ 0          │ J_PR,groR  │ J_PR,maiR  │ 0          │ 0          │ 0          ┃
┃VS  ┃ 0          │ J_VR,groR  │ 0          │ 0          │ 0          │ 0          ┃
┃RS  ┃ 0          │ 0          │ 0          │ 0          │ 0          │ 0          ┃
┃ECS ┃ J_ECR,assR │ J_ECR,groR │ J_ECR,maiR │ 0          │ J_ECS,rejR │ J_ECR,traR ┃
┃ENS ┃ J_ENR,assR │ J_ENR,groR │ J_ENR,maiR │ 0          │ J_ENS,rejR │ J_ENR,traR ┃
┃ES  ┃ J_ER,assR  │ J_ER,groR  │ J_ER,maiR  │ 0          │ 0          │ J_ER,traR  ┃
┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

J1: Catabolic flux diagram.
┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃    ┃ catS       │ rejS      │ losS mols!┃    ┃ catR       │ rejR      │ losR mols!┃
┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃EES ┃ J_EES,catS │ 0         │ 0         ┃EES ┃ J_EER,catR │ 0         │ 0         ┃
┃CNS ┃ J_CNS,catS │ 0         │ 0         ┃NS  ┃ J_CNR,catR │ 0         │ 0         ┃
┃CS  ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃CS  ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
┃NS  ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃ENS ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
┃ES  ┃ J_ES,catS  │ 0         │ 0         ┃ES  ┃ J_ER,catR  │ 0         │ 0         ┃
┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

M: State vector diagram.
┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
=#
