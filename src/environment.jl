export load_nichemap, apply_nichemap!, apply_deb_environment!, apply_energy_balance!

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
    # TODO: Separate out the photosynthetic spectrum.
    # env[:metout][:SOLR] = watts_to_light_mol.(env[:metout][:SOLR])

    return env
end

function apply_deb_environment!(settings, t::Real)
    update(s) = begin
        p = s.params
        p.J_L_F = radiation 
        p.X_H = soilmoist 
    end
    (airtemp, soiltemp, soilmoist, radiation, windspeed, humidity) = 
    get_nichemap(settings.environment, t)
    apply(update, structures)
    temps = celcius_to_kelvin.((airtemp, soiltemp))
    map(correct_temps, s, temps)
end

function apply_energy_balance!(settings, t::Real)
    update(s, airtemp, soiltemp, soilmoist, radiation, windspeed, humidity) = begin
      p = s.params
      p.photo_flux_density = radiation  
      p.Tair = airtemp
      p.rel_humidity = humidity
      p.wind = windspeed
      p.X_H = soilmoist 
    end
    structures = settings.structures
    updates = get_nichemap(settings.environment, t)
    # Update params for each structure
    apply(update, structures, updates...)

    (airtemp, soiltemp, soilmoist, radiation, windspeed, humidity) = updates
    setval!(structures[1].assim_state, structures[1].params)
    gas_ex!(structures[1].assim_state, structures[1].params)
    leaftemp = (structures[1].assim_state.Tleaf)
    temps = celcius_to_kelvin.((leaftemp, soiltemp))
    map(correct_temps!, structures, temps)
end

function get_nichemap(env::E, t::Real)::NTuple{6,Float64} where E
    t0 = 1 #/ settings.timestep_days
    pos = t0 + t
    airtemp = interpolate(env[:metout][:TALOC], pos) 
    soiltemp = interpolate(env[:soil][:D10cm], pos) 
    soilmoist = interpolate(env[:soilmoist][:WC10cm], pos)
    radiation = interpolate(env[:metout][:SOLR], pos) 
    windspeed = interpolate(env[:metout][:VLOC], pos) 
    humidity = interpolate(env[:humid][:RH20cm], pos) 

    return (airtemp, soiltemp, soilmoist, radiation, windspeed, humidity)
end

function correct_temps!(s::S, temp::Real) where S
    p = s.params
    corr = tempcorr(temp, p.REFERENCE_TEMP, p.ARRH_TEMP, p.LOWER_BOUNDARY, 
                    p.ARRH_LOWER, p.UPPER_BOUNDARY, p.ARRH_UPPER)
    # Scale params by temperature 
    for id in s.flags.temp
        s.params[id] = s.init_params[id] * corr
    end
end

interpolate(array::Array{Float64,1}, pos::Real)::Float64 = begin
    int = floor(Int64, pos)
    frac::Float64 = pos - int
    array[int] * (1.0 - frac) + array[int + 1] * frac
end
interpolate(array::Array{Float64,1}, pos::Int)::Float64 = array[int]
