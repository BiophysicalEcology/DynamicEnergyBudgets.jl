const parconv = 4.57mol*W^-1*s^-1

apply_environment!(plant, organs, u, t) = 
    apply_environment!(plant, environment(plant), organs, u, t)
apply_environment!(plant, ::Nothing, organs, u, t) = nothing
apply_environment!(plant, env, organs, u, t) = begin
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

apply_environment!(a::FvCBPhotosynthesis, o, u, shootenv, rootenv) = begin
    va = assimilation_vars(o)

    va.tair = temp = airtemperature(shootenv)
    va.windspeed = windspeed(shootenv)
    va.rh = relhumidity(shootenv)
    va.rnet = radiation(shootenv)
    va.par = radiation(shootenv) * parconv
    # va.soilmoist = mean_soilwatercontent(rootenv)
    va.swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)

    update_temp!(o, temp)

    if is_germinated(o, u)
        run_enbal!(assimilation_pars(o).photoparams, va)
    else
        va.tleaf = va.tair
    end
end

apply_environment!(a::AbstractCAssim, o, u, shootenv, rootenv) = begin
    update_temp!(o, airtemperature(shootenv))
    assimilation_vars(o).J_L_F = radiation(shootenv) * parconv
    # Use the mean soil water potential between zero and the full root depth 
    assimilation_vars(o).soilwaterpotential = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
end

apply_environment!(a::AbstractNAssim, o, u, rootenv, shootenv) = begin
    update_temp!(o, soiltemperature(rootenv))
    # va.X_H = interp_layer(inc_interp, env.soilwatercontent, i)
end

update_temp!(o, temp) = begin
    set_temp!(o, temp)
    set_tempcorrection!(o, tempcorr(temp, tempcorr_pars(o)))
end
