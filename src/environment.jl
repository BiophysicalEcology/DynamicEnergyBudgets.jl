const water_fraction_to_M = 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"

get_environment(::Type{Val{:soiltemperature}}, env::MicroclimateData, h, i) = 
    layer_interpolate(interp, env.soil, i) * u"°C"
get_environment(::Type{Val{:watercontent}}, env::MicroclimateData, h, i) =
    layer_interpolate(interp, env.soilmoist, i)
get_environment(::Type{Val{:waterpotential}}, env::MicroclimateData, h, i) = 
    layer_interpolate(interp, env.soilmoist, i) * u"Pa"
get_environment(::Type{Val{:airtemperature}}, env::MicroclimateData, h, i) =
    lin_interpolate(env.metout[:TALOC], i) * u"°C"
get_environment(::Type{Val{:windspeed}}, env::MicroclimateData, h, i) =
    lin_interpolate(env.metout[:VLOC], i) * u"m*s^-1"
get_environment(::Type{Val{:relhumidity}}, env::MicroclimateData, h, i) =
    layer_interpolate(interp, env.humid, i)
get_environment(::Type{Val{:radition}}, env::MicroclimateData, h, i) =
    lin_interpolate(env.metout[:SOLR], i) * u"W*m^-2"
get_environment(::Type{Val{:par}}, env::MicroclimateData, h, i) =
    lin_interpolate(env.metout[:SOLR], i) * 4.57u"mol*m^-2*s^-1"


apply_environment!(o, env, t) =
    apply_environment!(o, o.params.assimilation, env, t)

apply_environment!(o, a::AbstractAssimilation, env::Void, t) = nothing

apply_environment!(o, a::AbstractCarbonAssimilation, env, t) = begin
    p, v, u = components(o); va = v.assimilation;
    pos = ustrip(t) + 1
    h = v.height = allometric_height(p.allometry, o)
    interp = layer_setup(o.vars.height)

    va.tair = get_environment(Val{:airtemperature}, env, h, pos)
    va.windspeed = get_environment(Val{:windspeed}, env, h, pos)
    va.rh = get_environment(Val{:relhumidity}, env, h, pos)
    va.rnet = get_environment(Val{:radiation}, env, h, pos)
    va.par = get_environment(Val{:par}, env, h, pos)
    va.soilmoist = get_environment(Val{:soilwatercontent}, env, h, pos)
    va.swp = get_environment(Val{:soilwaterpotential}, env, h, pos)

    if germinated(u.V, p.M_Vgerm)
        phototranspiration!(va, va.photoparams)
    else
        va.tleaf = va.tair
    end

    correct_temps!(o, va.tleaf)
end

apply_environment!(o, a::KooijmanSLAPhotosynthesis, env, t) = begin
    p, v, u = components(o); va = v.assimilation;
    pos = ustrip(t) + 1
    h = v.height = allometric_height(p.allometry, o)

    va.tair = get_environment(Val{:airtemperature}, env, h, pos)
    va.J_L_F = get_environment(Val{:par}, env, h, pos)
    correct_temps!(o, va.tair)
end

apply_environment!(o, a::AbstractNitrogenAssimilation, env, t) = begin
    p, v, u = components(o); va = v.assimilation;
    pos = ustrip(t) + 1
    h = v.height = allometric_height(p.allometry, o)
    interp = layer_setup(v.height)
    va.temp = get_environment(Val{:soiltemperature}, env, h, pos)
    va.X_H = get_environment(Val{:soilwatercontent}, env, h, pos) * water_fraction_to_M

    correct_temps!(o, va.temp)
end

"Scale variables by temperature"
function correct_temps!(o, temp)
    p, v, u, sh = components(o)
    corr = tempcorr(temp, sh.tempcorr)
    v.k_E = p.k_E * corr
    v.k_EC = p.k_EC * corr
    v.k_EN = p.k_EN * corr
    v.j_E_mai = p.j_E_mai * corr
    v.j_E_rep_mai = p.maturity.j_E_rep_mai * corr
    v.j_P_mai = p.j_P_mai * corr
end
