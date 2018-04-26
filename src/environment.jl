
function load_nichemap(location, years::Real)
    # TODO: also get soil nitrate/ammonia levels.
    # TODO: Use NicheMapRs plant wilting output.

    # location = "Brisbane, qld"
    # years = 1
    env = nichemap_global(location, years)

    # Convert units
    for key in keys(env[:soilmoist])
        env[:soilmoist][key] = water_content_to_mols_per_litre.(env[:soilmoist][key])
    end

    return env
end

function apply_environment!(organ, 
                            env::NicheMapMicroclimate, 
                            t::Number)
    a = organ.values.assim_state
    pos = t + 1

    a.windspeed = interpolate(env.metout[:VLOC], pos) * u"m*s^-1"
    a.rh = interpolate(env.humid[:RH20cm], pos) 
    a.tair = interpolate(env.metout[:TALOC], pos) * u"°C"
    radiation = interpolate(env.metout[:SOLR], pos) 
    a.rnet = radiation * u"W*m^-2"
    a.par = radiation * 4.57u"mol*m^-2*s^-1" 

    photosynthesis_transpiration!(a, organ.params.assimilation.photoparams)
    correct_temps!(organ, a.tleaf)
end

function apply_environment!(organ::Organ, 
                            environment::NicheMapMicroclimate,
                            t::Number)
    a = organ.values.assim_state
    pos = t + 1

    a.soiltemp = interpolate(env.soil[:D10cm], pos) * u"°C"
    a.soilmoist = interpolate(env.soilmoist[:WC10cm], pos)
    a.soilpot = interpolate(env.soilmoist[:WC10cm], pos)
    a.soilpotshade = interpolate(env.shadpot[:PT10cm], pos)
    a.soilpotshade = interpolate(env.soilshade[:PT10cm], pos)

    correct_temps!(organ, a.soiltemp)
end

# function apply_environment!(environment::NicheMapMicroclimate, organ::Organ, t::Number)
#      update(s) = begin
#          p = s.params
#          p.J_L_F = radiation 
#          p.X_H = soilmoist 
#      end
#      (airtemp, soiltemp, soilmoist, radiation, windspeed, humidity) = 
#          get_nichemap(environment, t)
#      apply(update, organs)
#      temps = celcius_to_kelvin.((airtemp, soiltemp))
#      map(correct_temps, settings.structures, temps)
# end

function correct_temps!(o, temp::Number)
    p = o.params
    corr = tempcorr(temp |> u"K", getfield.(tc, fieldnames(tc))...)
    # Scale params by temperature 
    # for key in o.flags.temp
        value = getfield(o.init_params, key) * corr
        setfield!(o, key, value)
    # end
end

interpolate(array, pos::Number) = begin
    int = floor(Int64, pos)
    frac::Float64 = pos - int
    array[int] * (1 - frac) + array[int + 1] * frac
end
interpolate(array, pos::Int) = array[int]
