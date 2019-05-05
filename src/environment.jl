const parconv = 4.57mol*W^-1*s^-1

struct ManualTemperature{T}
    airtemperature::T
    soiltemperature::T
end

apply_environment!(plant, organs, u, t) = 
    apply_environment!(plant, environment(plant), organs, u, t)
apply_environment!(plant, env::Nothing, organs, u, t) = nothing
apply_environment!(plant, env::ManualTemperature, organs, u, t) = begin
    envpos = ustrip(calc_envtime(plant, t))
    shoot, root = organs
    update_temp!(shoot, env.airtemperature)
    update_temp!(root, env.soiltemperature)
end

apply_environment!(plant, env::AbstractMicroclimate, organs, u, t) = begin
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

apply_environment!(plant, env::MicroclimControl, organs, u, t) = begin
    shoot, root = organs
    apply_environment!(assimilation_pars(shoot), shoot, u, env, env)
    apply_environment!(assimilation_pars(root), root, u, env, env)
    nothing
end

apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv) = begin
    v = a.vars

    v.tair = temp = airtemperature(shootenv)
    v.windspeed = windspeed(shootenv)
    v.rh = relhumidity(shootenv)
    v.rnet = radiation(shootenv)
    v.par = radiation(shootenv) * parconv
    v.vpd = vapour_pressure_deficit(v.tair, v.rh)
    # v.soilmoist = mean_soilwatercontent(rootenv)
    v.swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)

    enbal!(assimilation_pars(o).photoparams, v)

    update_temp!(o, v.tleaf)
end

apply_environment!(a::AbstractAssim, o, u, shootenv, rootenv) = 
    update_temp!(o, airtemperature(shootenv))
apply_environment!(a::AbstractNAssim, o, u, rootenv, shootenv) = 
    update_temp!(o, soiltemperature(rootenv))

apply_environment!(a::AbstractKooijmanPhoto, o, u, shootenv, rootenv) = begin
    update_temp!(o, airtemperature(shootenv))
    a.vars.J_L_F = radiation(shootenv) * parconv
    # Use the mean soil water potential between zero and the full root depth 
    a.vars.soilwaterpotential = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
end


update_temp!(o, temp) = begin
    set_temp!(o, temp)
    set_tempcorrection!(o, tempcorr(tempcorr_pars(o), temp))
end


