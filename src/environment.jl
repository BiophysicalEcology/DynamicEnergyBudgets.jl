
const nichemap_increments = (0.0u"cm",2.5u"cm",5.0u"cm",10.0u"cm",15.0u"cm", 20.0u"cm",30.0u"cm",50.0u"cm",100.0u"cm",200.0u"cm")

struct Interp
    lower::Int
    upper::Int
    lowerfrac::Float64
    upperfrac::Float64
end

allometric_height(f::SqrtAllometry, p, u) = 
    sqrt((u.P * p.w_P + u.V * p.w_V) / oneunit(u.V*p.w_V)) * f.scale

height_interpolate(height) = begin
    for (i, upper) in enumerate(nichemap_increments)
        if upper > height
            lower = nichemap_increments[i - 1]
            p = (height-lower)/(upper-lower)
            return Interp(i+1, i+2, p, 1.0-p)
        end
    end
    # Otherwise its taller/deeper than we have data, use the largest we have.
    max = length(nichemap_increments) + 2
    return Interp(max, max, 1.0, 0.0)
end

apply_environment!(organ::Organ, env, t) =
    apply_environment!(organ, organ.params.assimilation, env, t)

function apply_environment!(o::Organ, a::CarbonAssimilation, env::NicheMapMicroclimate, t)
    v = o.vars.assimilation
    pos = ustrip(t) + 1
    o.vars.height = allometric_height(o.params.allometry, o.params, o.state)
    interp = height_interpolate(o.vars.height)

    v.windspeed = interpolate(env.metout[:VLOC], pos) * u"m*s^-1"
    v.rh = combine(interp, env.humid, pos)
    v.tair = interpolate(env.metout[:TALOC], pos) * u"°C"
    radiation = interpolate(env.metout[:SOLR], pos)
    v.rnet = radiation * u"W*m^-2"
    v.par = radiation * 4.57u"mol*m^-2*s^-1"
    v.soilmoist = combine(interp, env.soilmoist, pos)
    v.swp = combine(interp, env.soilpot, pos) * u"kPa"
    phototranspiration!(v, a.photoparams)
    correct_temps!(o, v.tleaf)
end

function apply_environment!(o::Organ, a::NitrogenAssimilation, env::NicheMapMicroclimate, t)
    v = o.vars.assimilation
    pos = ustrip(t) + 1
    o.vars.height = allometric_height(o.params.allometry, o.params, o.state)
    interp = height_interpolate(o.vars.height)

    v.temp = combine(interp, env.soil, pos) * u"°C"
    v.soilmoist = combine(interp, env.soilmoist, pos)
    v.soilpot = combine(interp, env.soilmoist, pos)
    v.soilpotshade = combine(interp, env.shadpot, pos)
    correct_temps!(o, v.temp)
end

"Scale variables by temperature"
function correct_temps!(o, temp)
    v = o.vars; p = o.params
    corr = tempcorr(temp, p.tempcorr)
    v.k_E = p.k_E * corr
    v.k_EC = p.k_EC * corr
    v.k_EN = p.k_EN * corr
    v.j_E_mai = p.j_E_mai * corr
    v.j_E_rep_mai = p.j_E_rep_mai * corr
    v.j_P_mai = p.j_P_mai * corr
end

combine(i, data, pos) = begin
    interpolate(data[i.lower], pos) * i.lowerfrac +
    interpolate(data[i.upper], pos) * i.upperfrac
end

interpolate(array, pos::Number) = begin
    int = floor(Int64, pos)
    frac::Float64 = pos - int
    array[int] * (1 - frac) + array[int + 1] * frac
end
interpolate(array, pos::Int) = array[pos]
