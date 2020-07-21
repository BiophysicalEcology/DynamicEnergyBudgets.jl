using .Microclimate

"""
    apply_environment!(plant::Plant, env::AbstractMicroclimate, organs, u, t)

Apply a microclimate from Microclimates.jl to a `Plant`.
"""
apply_environment!(env::AbstractMicroclimate, plant::Plant, organs, u, t) = begin
    envind = calc_envind(plant, t)
    shoot, root = organs
    # Get the climate at a particular height and time
    # Use half the height for the median temperature/windspeed/humidity etc.
    # relating temperature to the distribution of biomass may be more accurate.
    shootenv = MicroclimInstant(env, envind, height(shoot)/2)
    rootenv = MicroclimInstant(env, envind, depth(root)/2)
    # Update any environment dependent variables
    apply_environment!(assimilation_pars(shoot), shoot, u, shootenv, rootenv)
    apply_environment!(assimilation_pars(root), root, u, rootenv, shootenv)
    nothing
end

# We need an environment data point for t = 0 even if we
# never use it - we may need to interpolate at t = 0.5
calc_envind(o, t) = ustrip(o.environment_start[] + t) + 1

"""
    apply_environment!(plant::Plant, env::MicroclimControl, organs, u, t)

Apply an environment to a Plant using live user input for variables,
instead of data.
"""
apply_environment!(env::MicroclimControl, plant::Plant, organs, u, t) = begin
    shoot, root = organs
    apply_environment!(assimilation_pars(shoot), shoot, u, env, env)
    apply_environment!(assimilation_pars(root), root, u, env, env)
    nothing
end


"""
    apply_environment!(a::AbstractKooijmanPhoto, o, u, shootenv, rootenv)

Apply environment to shoot with `AbstractKooijmanPhoto` assimilation model.
"""
apply_environment!(a::AbstractKooijmanPhoto, o, u, shootenv, rootenv) = begin
    v = a.vars
    update_temp!(o, airtemperature(shootenv))
    v.J_L_F = radiation(shootenv) * parconv
    # Use the mean soil water potential between zero and the full root depth 
    # swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
    swp = if typeof(rootenv) <: MicroclimControl
        soilwaterpotential(rootenv)
    else
        layermax(soilwaterpotential(rootenv.microclimate), rootenv)
    end
    v.soilwaterpotential = swp
    set_swp!(o, swp)
end
"""

Apply environment to shoot with any `AbstractAssim` assimilation model.
This simply adds the air temperature.
"""
apply_environment!(a::AbstractAssim, o, u, shootenv, rootenv) = 
    update_temp!(o, airtemperature(shootenv))
"""
    apply_environment!(a::AbstractNAssim, o, u, rootenv, shootenv)

Apply environment to root with `AbstractNAssim` assimilation model.
This simply adds the soil temperature.
"""
apply_environment!(a::AbstractNAssim, o, u, rootenv, shootenv) = 
    update_temp!(o, soiltemperature(rootenv))

