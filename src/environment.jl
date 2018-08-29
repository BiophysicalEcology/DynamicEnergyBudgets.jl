const water_fraction_to_M = 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"

get_environment(t::Type{Val{:soiltemperature}}, env::M, interp, i) where M <: MicroclimateTable = 
    layer_interp(interp, env.soil, i) * u"°C"
get_environment(t::Type{Val{:soilwatercontent}}, env::M, interp, i) where M <: MicroclimateTable =
    layer_interp(interp, env.soilmoist, i)
get_environment(t::Type{Val{:soilwaterpotential}}, env::M, interp, i) where M <: MicroclimateTable = 
    layer_interp(interp, env.soilmoist, i) * u"Pa"
get_environment(t::Type{Val{:relhumidity}}, env::M, interp, i) where M <: MicroclimateTable =
    layer_interp(interp, env.humid, i)

get_environment(t::Type{Val{:airtemperature}}, env::M, interp, i) where M <: MicroclimateTable =
    lin_interp(env.metout, Val{:TALOC}, i) * u"°C"
get_environment(t::Type{Val{:windspeed}}, env::M, interp, i) where M <: MicroclimateTable =
    lin_interp(env.metout, Val{:VLOC}, i) * u"m*s^-1"
get_environment(t::Type{Val{:radiation}}, env::M, interp, i) where M <: MicroclimateTable =
    lin_interp(env.metout, Val{:SOLR}, i) * u"W*m^-2"
get_environment(t::Type{Val{:par}}, env::M, interp, i) where M <: MicroclimateTable =
    lin_interp(env.metout, Val{:SOLR}, i) * 4.57u"mol*m^-2*s^-1"


apply_environment!(organs::Tuple{O,Vararg}, ux::Tuple{U,Vararg}, env, t) where {U,O} = begin
    env == nothing && return 
    apply_environment!(organs[1], ux[1], env, t)
    apply_environment!(tail(organs), tail(ux), env, t)
    nothing
end
apply_environment!(organs::Tuple{}, ux::Tuple{}, env, t) = nothing

apply_environment!(o::Organ, u, env, t) = begin
    apply_environment!(o.params.assimilation, o, u, env, t)
    nothing
end
apply_environment!(o::Organ, u, env::Nothing, t) = nothing

apply_environment!(a::FvCBPhotosynthesis, o, u, env, t) = begin
    p, v = unpack(o); va = assimilation(v);
    pos = ustrip(t) + 1
    h = allometric_height(p.allometry, o, u)
    set_var!(v, :height, h)
    interp = LayerInterp(height(v))

    va.tair = get_environment(Val{:airtemperature}, env, interp, pos)
    va.windspeed = get_environment(Val{:windspeed}, env, interp, pos)
    va.rh = get_environment(Val{:relhumidity}, env, interp, pos)
    va.rnet = get_environment(Val{:radiation}, env, interp, pos)
    va.par = get_environment(Val{:par}, env, interp, pos)
    va.soilmoist = get_environment(Val{:soilwatercontent}, env, interp, pos)
    va.swp = get_environment(Val{:soilwaterpotential}, env, interp, pos)

    if germinated(u[V], p.M_Vgerm)
        phototranspiration!(va, p.assimilation.photoparams)
    else
        va.tleaf = temp(v)
    end

    set_var!(v, :temp, va.tleaf)
    set_var!(v, :tempcorrection, tempcorr(temp(v), o.shared.tempcorr))
    nothing
end

apply_environment!(a::AbstractCAssim, o, u, env, t) = begin
    p, v = unpack(o); va = assimilation(v);
    pos = ustrip(t) + 1
    h = allometric_height(p.allometry, o, u)
    set_var!(v, :height, h)
    interp = LayerInterp(h)

    set_var!(v, :temp, get_environment(Val{:airtemperature}, env, interp, pos))
    va.J_L_F = get_environment(Val{:par}, env, interp, pos)
    set_var!(v, :tempcorrection,  tempcorr(temp(v), o.shared.tempcorr))
    nothing
end

apply_environment!(a::Nothing, o, u, env, t) = begin
    p, v = unpack(o); va = assimilation(v);
    pos = ustrip(t) + 1
    h = allometric_height(p.allometry, o, u)
    set_var!(v, :height, h)
    interp = LayerInterp(h)

    set_var!(v, :temp, get_environment(Val{:airtemperature}, env, interp, pos))
    set_var!(v, :tempcorrection,  tempcorr(temp(v), o.shared.tempcorr))
    nothing
end

apply_environment!(a::AbstractNAssim, o, u, env, t) = begin
    p, v = unpack(o); va = assimilation(v);
    pos = ustrip(t) + 1
    h = allometric_height(p.allometry, o, u)
    set_var!(v, :height, h)
    interp = LayerInterp(h)

    set_var!(v, :temp, get_environment(Val{:soiltemperature}, env, interp, pos))
    va.X_H = get_environment(Val{:soilwatercontent}, env, interp, pos) * water_fraction_to_M

    set_var!(v, :tempcorrection,  tempcorr(temp(v), o.shared.tempcorr))
    nothing
end
