const parconv = 4.57mol*W^-1*s^-1

"""
    ManualTemperature(airtemperature, soiltemperature)

Environment for simple manual temperature control.
"""
struct ManualTemperature{T}
    airtemperature::T
    soiltemperature::T
end

"""
    apply_environment!(organism::Organism, organs::Tuple, u, t) 

Apply environmental conditions to an organism. 
"""
apply_environment!(organism::AbstractOrganism, organs, u, t) = 
    apply_environment!(organism, environment(organism), organs, u, t)

"""
    apply_environment!(organism, env::Nothing, organs, u, t)

No environment included, do nothing. Initial environment variables
will be used for the whole simulation
"""
apply_environment!(organism, env::Nothing, organs, u, t) = nothing

"""
    apply_environment!(plant, env::ManualTemperature, organs, u, t)

Apply an environment to a plant using live user input for variables,
instead of data.
"""
apply_environment!(plant, env::ManualTemperature, organs, u, t) = begin
    envpos = ustrip(calc_envtime(plant, t))
    shoot, root = organs
    update_temp!(shoot, env.airtemperature)
    update_temp!(root, env.soiltemperature)
end

"""
    apply_environment!(plant::Plant, env::AbstractMicroclimate, organs, u, t)

Apply a microclimate from Microclimates.jl to a `Plant`.
"""
apply_environment!(plant::Plant, env::AbstractMicroclimate, organs, u, t) = begin
    envpos = ustrip(calc_envtime(plant, t))
    shoot, root = organs
    # Get the climate at a particular height and time
    # Use half the height for the median temperature/windspeed/humidity etc.
    # relating temperature to the distribution of biomass may be more accurate.
    shootenv = MicroclimInstant(env, height(root)/2, envpos)
    rootenv = MicroclimInstant(env, depth(root)/2, envpos)
    # Update any environment dependent variables
    apply_environment!(assimilation_pars(shoot), shoot, u, shootenv, rootenv)
    apply_environment!(assimilation_pars(root), root, u, rootenv, shootenv)
    nothing
end

"""
    apply_environment!(plant::Plant, env::MicroclimControl, organs, u, t)

Apply an environment to a Plant using live user input for variables,
instead of data.
"""
apply_environment!(plant::Plant, env::MicroclimControl, organs, u, t) = begin
    shoot, root = organs
    apply_environment!(assimilation_pars(shoot), shoot, u, env, env)
    apply_environment!(assimilation_pars(root), root, u, env, env)
    nothing
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

"""
    apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv)

Apply environment to shoot with `AbstractFvCBCAssim` assimilation model.
"""
apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv) = begin
    v = a.vars

    v.tair = airtemperature(shootenv)
    v.windspeed = windspeed(shootenv)
    v.rh = relhumidity(shootenv)
    v.rnet = radiation(shootenv)
    v.par = radiation(shootenv) * parconv
    v.vpd = vapour_pressure_deficit(v.tair, v.rh)
    # v.soilmoist = mean_soilwatercontent(rootenv)
    # swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
    swp = if typeof(rootenv) <: MicroclimControl
        soilwaterpotential(rootenv)
    else
        layermax(soilwaterpotential(rootenv.microclimate), rootenv)
    end
    v.swp = swp
    set_swp!(o, v.swp)

    # This runs energy balance which contains photosynthesis
    # Later in assimilation we can just use the result
    enbal!(v, assimilation_pars(o).photoparams)
    set_soilcorrection!(o, v.fsoil)

    update_temp!(o, v.tleaf)
end

"""
    apply_environment!(a::AbstractKooijmanPhoto, o, u, shootenv, rootenv)

Apply environment to shoot with `AbstractKooijmanPhoto` assimilation model.
"""
apply_environment!(a::AbstractKooijmanPhoto, o, u, shootenv, rootenv) = begin
    update_temp!(o, airtemperature(shootenv))
    a.vars.J_L_F = radiation(shootenv) * parconv
    # Use the mean soil water potential between zero and the full root depth 
    # swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
    swp = layermax(soilwaterpotential(rootenv.microclimate), rootenv)
    a.vars.soilwaterpotential = swp
    set_swp!(o, swp)
end

"""
    update_temp!(o, temp)

Set the temperature, and update the temperature correction variable.
"""
update_temp!(o, temp) = begin
    set_temp!(o, temp)
    set_tempcorrection!(o, tempcorr(tempcorr_pars(o), temp))
end


