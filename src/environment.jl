const water_fraction_to_M = 1.0m^3*m^-3 * 1kg*L^-1 / 18.0g*mol^-1

get_environment(t::Type{Val{:soiltemperature}}, env::M, interp, i) where M <: MicroclimateTables = 
    layer_interp(interp, env.soil, i) * °C
get_environment(t::Type{Val{:soilwatercontent}}, env::M, interp, i) where M <: MicroclimateTables =
    layer_interp(interp, env.soilmoist, i)
get_environment(t::Type{Val{:soilwaterpotential}}, env::M, interp, i) where M <: MicroclimateTables = 
    layer_interp(interp, env.soilmoist, i) * Pa
get_environment(t::Type{Val{:relhumidity}}, env::M, interp, i) where M <: MicroclimateTables =
    layer_interp(interp, env.humid, i)

get_environment(t::Type{Val{:airtemperature}}, env::M, interp, i) where M <: MicroclimateTables =
    lin_interp(env.metout.TALOC, i) * °C
get_environment(t::Type{Val{:windspeed}}, env::M, interp, i) where M <: MicroclimateTables =
    lin_interp(env.metout.VLOC, i) * m*s^-1
get_environment(t::Type{Val{:radiation}}, env::M, interp, i) where M <: MicroclimateTables =
    lin_interp(env.metout.SOLR, i) * W*m^-2
get_environment(t::Type{Val{:par}}, env::M, interp, i) where M <: MicroclimateTables =
    lin_interp(env.metout.SOLR, i) * 4.57mol*m^-2*s^-1


apply_environment!(organs::Tuple{O,Vararg}, ux::Tuple{U,Vararg}, env, t) where {U,O} = begin
    env == nothing && return 
    apply_environment!(organs[1], ux[1], env, t)
    apply_environment!(tail(organs), tail(ux), env, t)
end
apply_environment!(organs::Tuple{}, ux::Tuple{}, env, t) = nothing

apply_environment!(o::Organ, u, env::Nothing, t) = nothing
apply_environment!(o::Organ, u, env, t) = begin
    pos = ustrip(t) + 1
    interp = LayerInterp(height(o.vars))
    apply_environment!(assimilation_pars(o), o, u, env, t, pos, interp)
    set_var!(o, :tempcorrection, tempcorr(o))
    nothing
end

apply_environment!(a::FvCBPhotosynthesis, o, u, env, t, pos, interp) = begin
    va = assimilation_vars(o.vars)

    va.tair = get_environment(Val{:airtemperature}, env, interp, pos)
    va.windspeed = get_environment(Val{:windspeed}, env, interp, pos)
    va.rh = get_environment(Val{:relhumidity}, env, interp, pos)
    va.rnet = get_environment(Val{:radiation}, env, interp, pos)
    va.par = get_environment(Val{:par}, env, interp, pos)
    va.soilmoist = get_environment(Val{:soilwatercontent}, env, interp, pos)
    va.swp = get_environment(Val{:soilwaterpotential}, env, interp, pos)

    if germinated(u, M_Vgerm(o))
        phototranspiration!(va, assimilation_pars(o).photoparams)
    else
        va.tleaf = temp(v)
    end
end

apply_environment!(a::AbstractCAssim, o, u, env, t, pos, interp) = begin
    set_var!(o, :temp, get_environment(Val{:airtemperature}, env, interp, pos))
    assimilation_vars(o.vars).J_L_F = get_environment(Val{:par}, env, interp, pos)
end

apply_environment!(a::AbstractNAssim, o, u, env, t, pos, interp) = begin
    set_var!(o, :temp, get_environment(Val{:soiltemperature}, env, interp, pos))
    # va.X_H = get_environment(Val{:soilwatercontent}, env, interp, pos) * water_fraction_to_M
end
