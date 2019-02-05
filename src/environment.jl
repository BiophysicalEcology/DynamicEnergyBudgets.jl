const parconv = 4.57mol*W^-1*s^-1

apply_environment!(plant, organs, u, t) = 
    apply_environment!(plant, environment(plant), organs, u, t)
apply_environment!(plant, ::Nothing, organs, u, t) = nothing
apply_environment!(plant, env, organs, u, t) = begin
    envpos = ustrip(calc_envtime(plant, t))
    shoot, root = organs
    # Get the climate at a particular height and time
    shoot_env = MicroclimInstant(env, height(shoot), envpos)
    root_env = MicroclimInstant(env, height(root), envpos)
    # Update any environment dependent variables
    apply_environment!(assimilation_pars(shoot), shoot, u, shoot_env, root_env)
    apply_environment!(assimilation_pars(root), root, u, root_env, shoot_env)
    nothing
end

apply_environment!(a::FvCBPhotosynthesis, o, u, shootenv, rootenv) = begin
    va = assimilation_vars(o)

    va.tair = temp = airtemperature(env)
    va.windspeed = windspeed(env)
    va.rh = relhumidity(env)
    va.rnet = radiation(env)
    va.par = radiation(env) * parconv
    va.soilmoist = soilwatercontent(env)
    va.swp = mean_soilwaterpotential(env)

    update_temp!(o, temp)

    if is_germinated(o, u)
        phototranspiration!(va, assimilation_pars(o).photoparams)
    else
        va.tleaf = va.tair
    end
end

apply_environment!(a::AbstractCAssim, o, u, shootenv, rootenv) = begin
    temp = airtemperature(shootenv)
    update_temp!(o, temp)
    assimilation_vars(o).J_L_F = radiation(shootenv) * 4.57mol*W^-1*s^-1
    assimilation_vars(o).soilwaterpotential = mean_soilwaterpotential(rootenv)
end

apply_environment!(a::AbstractNAssim, o, u, rootenv, shootenv) = begin
    temp = soiltemperature(rootenv)
    update_temp!(o, temp)
    # va.X_H = interp_layer(inc_interp, env.soilwatercontent, i)
end

update_temp!(o, temp) = begin
    set_temp!(o, temp)
    set_tempcorrection!(o, tempcorr(temp, tempcorr_pars(o)))
end
