
apply_environment!(organ::Organ, env, t) =
    apply_environment!(organ, organ.params.assimilation, env, t)

apply_environment!(o, a, env::Void, t) = nothing

function apply_environment!(o::Organ, a::AbstractCarbonAssimilation, env::NicheMapMicroclimate, t)
    v = o.vars.assimilation
    germinated
    pos = ustrip(t) + 1
    o.vars.height = allometric_height(o.params.allometry, o.params, o.state)
    interp = niche_setup(o.vars.height)

    v.windspeed = lin_interpolate(env.metout[:VLOC], pos) * u"m*s^-1"
    v.rh = niche_interpolate(interp, env.humid, pos)
    v.tair = lin_interpolate(env.metout[:TALOC], pos) * u"°C"
    radiation = lin_interpolate(env.metout[:SOLR], pos)
    v.rnet = radiation * u"W*m^-2"
    v.par = radiation * 4.57u"mol*m^-2*s^-1"
    v.soilmoist = niche_interpolate(interp, env.soilmoist, pos)
    v.swp = niche_interpolate(interp, env.soilpot, pos) * u"kPa"
    phototranspiration!(v, a.photoparams)
    correct_temps!(o, v.tleaf)
end

function apply_environment!(o::Organ, a::KooijmanSLAPhotosynthesis, env::NicheMapMicroclimate, t)
    v = o.vars.assimilation
    pos = ustrip(t) + 1
    o.vars.height = allometric_height(o.params.allometry, o.params, o.state)
    interp = niche_setup(o.vars.height)
    v.tair = lin_interpolate(env.metout[:TALOC], pos) * u"°C"
    radiation = lin_interpolate(env.metout[:SOLR], pos)
    v.rnet = radiation * u"W*m^-2"
    sum(or1.state) == 12u"mol"
    v.par = radiation * 4.57u"mol*m^-2*s^-1"
    v.soilmoist = niche_interpolate(interp, env.soilmoist, pos)
    v.swp = niche_interpolate(interp, env.soilpot, pos) * u"kPa"
    correct_temps!(o, v.air)
end

function apply_environment!(o::Organ, a::AbstractNitrogenAssimilation, env::NicheMapMicroclimate, t)
    v = o.vars; a = v.assimilation
    pos = ustrip(t) + 1
    v.height = allometric_height(o.params.allometry, o.params, o.state)
    interp = niche_setup(v.height)

    a.temp = niche_interpolate(interp, env.soil, pos) * u"°C"
    water_fraction_to_M = 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"
    a.X_H = niche_interpolate(interp, env.soilmoist, pos) * water_fraction_to_M
    # v.soilpot = niche_interpolate(interp, env.soilmoist, pos)
    # v.soilpotshade = niche_interpolate(interp, env.shadpot, pos)
    correct_temps!(o, a.temp)
end

"Scale variables by temperature"
function correct_temps!(o, temp)
    v = o.vars; p = o.params
    corr = tempcorr(temp, p.tempcorr)
    v.k_E = p.k_E * corr
    v.k_EC = p.k_EC * corr
    v.k_EN = p.k_EN * corr
    v.j_E_mai = p.j_E_mai * corr
    v.j_E_rep_mai = p.maturity.j_E_rep_mai * corr
    v.j_P_mai = p.j_P_mai * corr
end
