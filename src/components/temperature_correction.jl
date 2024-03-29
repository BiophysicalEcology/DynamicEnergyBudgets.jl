"""
Temperature correction parameters
"""
abstract type AbstractTemperatureCorrection{T} end

# Temperature mixins
@mix @columns struct Tbase{T}
    # Field       | Default | Unit | Bounds          | Log | Description
    reftemp::T    | 300.0   | K    | (273.0, 325.0)  | _   | "Reference temperature for all rate parameters"
    arrtemp::T    | 2000.0  | K    | (200.0, 4000.0) | _   | "Arrhenius temperature"
end

@mix @columns struct Tlow{T}
    tbelow::T     | -30.0   | K | (0.0, -40.0)       | _ | "Lower boundary of tolerance range"
    arrlower::T   | 20000.0 | K | (2000.0, 40000.0)  | _ | "Arrhenius temperature for lower boundary"
end

@mix @columns struct Tup{T}
    tabove::T     | 5.0     | K | (0.0, 20)          | _ | "Upper boundary of tolerance range"
    arrupper::T   | 70000.0 | K | (7000.0, 140000.0) | _ | "Arrhenius temperature for upper boundary"
end

"""
    tempcorr(f, t)

Temperature related correction of growth rate.

-`f`: a formulation struct or `nothing` for no temperature correction.
-`t`: the current temperature, in degrees Celcius or Kelvin (Unitful.jl u"°C" or u"K")
"""
function tempcorr end
tempcorr(f::Nothing, t) = 1

"""
    TempCorr(reftemp, arrtemp)

Simple temperature correction parameters

$(FIELDDOCTABLE)
"""
@Tbase struct TempCorr{T} <: AbstractTemperatureCorrection{T} end

tempcorr(f::TempCorr, t) = tempcorr(t |> K, f.reftemp, f.arrtemp)

"""
    TempCorrLower(reftemp, arrtemp, tbelow, arrlower)

Temperature correction with lower bounds parameters.

$(FIELDDOCTABLE)
"""
@Tbase @Tlow struct TempCorrLower{T} <: AbstractTemperatureCorrection{T} end

tempcorr(f::TempCorrLower, t) =
    tempcorr(t |> K, f.reftemp, f.arrtemp, f.tbelow + f.reftemp, f.arrlower)

"""
    TempCorrLowerUpper(reftemp, arrtemp, tbelow, arrlower, tabove, arrupper)

Temperature correction with lower and upper bound parameters.

$(FIELDDOCTABLE)
"""
@Tbase @Tlow @Tup struct TempCorrLowerUpper{T} <: AbstractTemperatureCorrection{T} end

tempcorr(f::TempCorrLowerUpper, t) =
    tempcorr(t |> K, f.reftemp, f.arrtemp, f.tbelow + f.reftemp, f.arrlower, f.tabove + f.reftemp, f.arrupper)

tempcorr(t, tref, a) = exp(a/tref - a/t)
tempcorr(t, tref, a, l, al) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l)) / (1.0 + exp(al/t - al/l))
tempcorr(t, tref, a, l, al, h, ah) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l) + exp(ah/h - ah/tref)) /
    (1.0 + exp(al/t - al/l) + exp(ah/h - ah/t))

"""
    ParentTardieu(ΔH_A, α, t0)

Simple 3 parameter temperature correction method. 
Growth response to temperature has smoother transients in plants than in animals, 
and a simpler formulation is more applicable.

$(FIELDDOCTABLE)
"""
@flattenable @columns struct ParentTardieu{Δ,Al,T,A} <: AbstractTemperatureCorrection{T}
    ΔH_A::Δ | true  | 63.5  | kJ*mol^-1 | (55.0, 65.0)   | _ | "The enthalpy of activation of the reaction. Determines the curvature at low temperature"
    α::Al   | true  | 3.5   | _         | (1.0, 10.0)    | _ | "The ratio ΔH_D / ΔH_A"
    t0::T   | true  | 300.0 | K         | (273.0, 325.0) | _ | "Reference temperature"
    A::A    | false | 0.0   | K^-1      | (0.0, 1.0)     | _ | "Trait scaling coefficient, calculated from other params"
end          
# Allways recalculate A
ParentTardieu(ΔH_A, α, t0, A) = ParentTardieu(ΔH_A, α, t0)
# Precalculate scalar A using the reference temperature
# to normalise the temperature response with the reference at 1.0
ParentTardieu(ΔH_A, α, t0) = begin
    A = 1 / parent_tardieua_unscaled(ΔH_A, α, t0, t0)
    ParentTardieu{typeof.((ΔH_A, α, t0, A))...}(ΔH_A, α, t0, A)
end

parent_tardieua_unscaled(ΔH_A, α, t0, t) = begin
    ex = exp(-ΔH_A / (R * t))
    t * ex / (1 + ex^(α * (1 - (t / t0))))
end

tempcorr(f::ParentTardieu, t) = f.A * parent_tardieua_unscaled(f.ΔH_A, f.α, f.t0, t)

