abstract type AbstractAssimilation end
abstract type AbstractCarbonAssimilation <: AbstractAssimilation end
abstract type AbstractNitrogenAssimilation <: AbstractAssimilation end
abstract type AbstractNH4_NO3Assimilation <: AbstractNitrogenAssimilation end

@mix @label @range @with_kw struct SLA
    SLA::typeof(1.0u"m^2*g^-1") = 9.10u"m^2*g^-1" | [5.0, 30.0]u"m^2*g^-1" | "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"
end

@SLA mutable struct C3Photosynthesis{P} <: AbstractCarbonAssimilation
    photoparams::P = Photosynthesis.EnergyBalance() | _ | _
end

@SLA mutable struct KooijmanSLAPhotosynthesis <: AbstractCarbonAssimilation
    k_C_binding::typeof(1.0u"mol*mol^-1*s^-1") = 1.0u"mol*mol^-1*s^-1"                 | _ | "scaling rate for carbon dioxide"
    k_O_binding::typeof(1.0u"mol*mol^-1*s^-1") = 1.0u"mol*mol^-1*s^-1"                 | _ | "scaling rate for oxygen"
    K_C::typeof(1.0u"mol*L^-1")                = fraction_per_litre_gas_to_mols(40.0/1e6)u"mol*L^-1"  | _ | "half-saturation concentration of carbon dioxide"
    K_O::typeof(1.0u"mol*L^-1")                = fraction_per_litre_gas_to_mols(0.0021)u"mol*L^-1"    | _ | "half-saturation concentration of oxygen"
    J_L_K::typeof(1.0u"mol*m^-2*s^-1")         = watts_to_light_mol(300.0)u"mol*m^-2*s^-1" | _  | "half-saturation flux of useful photons"
    j_L_Amax::typeof(1.0u"μmol*m^-2*s^-1")     = 20.0u"μmol*m^-2*s^-1"                     | _  | "max spec uptake of useful photons"
    j_C_Amax::typeof(1.0u"μmol*m^-2*s^-1")     = 90.0u"μmol*m^-2*s^-1"                     | _  | "max spec uptake of carbon dioxide"
    j_O_Amax::typeof(1.0u"μmol*m^-2*s^-1")     = 0.001u"μmol*m^-2*s^-1"                    | _  | "max spec uptake of oxygen"
end

@label @range @with_kw mutable struct Kooijman_NH4_NO3Assimilation <: AbstractNH4_NO3Assimilation 
    j_NH_Amax::typeof(1.0u"μmol*mol^-1*s^-1")  = 50.0u"μmol*mol^-1*s^-1" | [0.1, 1000.0]u"μmol*mol^-1*s^-1" | "max spec uptake of ammonia"
    j_NO_Amax::typeof(1.0u"μmol*mol^-1*s^-1")  = 50.0u"μmol*mol^-1*s^-1" | [0.1, 1000.0]u"μmol*mol^-1*s^-1" | "max spec uptake of nitrate"
    ρNO::typeof(1.0)                           = 0.7                     | [0.0, 1.0]                       | "weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::typeof(1.0u"mol*mol^-1")        = 1.25u"mol*mol^-1"       | _                                | "from roots C-reserve to reserve using ammonia"
    K_NH::typeof(1.0u"mmol*L^-1")              = 10.0u"mmol*L^-1"        | _                                | "half-saturation concentration of ammonia"
    K_NO::typeof(1.0u"mmol*L^-1")              = 10.0u"mmol*L^-1"        | _                                | "half-saturation concentration of nitrate"
    K_H::typeof(1.0u"mol*L^-1")                = 10.0u"mol*L^-1"         | _                                | "half-saturation concentration of water"
end

@label @range @with_kw mutable struct N_Assimilation <: AbstractNitrogenAssimilation
    j_N_Amax::typeof(1.0u"μmol*mol^-1*s^-1")   = 50.0u"μmol*mol^-1*s^-1" | [0.1, 1000.0]u"μmol*mol^-1*s^-1" | "max spec uptake of ammonia"
    K_N::typeof(1.0u"mmol*L^-1")               = 10.0u"mmol*L^-1"        | _                                | "half-saturation concentration of nitrate"
    K_H::typeof(1.0u"mol*L^-1")                = 10.0u"mol*L^-1"         | _                                | "half-saturation concentration of water"
end


abstract type AbstractTempCorr end

@mix @label @range @with_kw struct Tbase
    reftemp::typeof(1.0u"K")    = 310.0u"K"   | [273.0, 325.0]u"K"     | "Reference temperature for all rate parameters"
    arrtemp::typeof(1.0u"K")    = 2000.0u"K"  | [200.0, 4000.0]u"K"    | "Arrhenius temperature"
end
@mix struct Tlow
    lowerbound::typeof(1.0u"K") = 280.0u"K"   | [273.0, 325.0]u"K"     | "Lower boundary of tolerance range"
    arrlower::typeof(1.0u"K")   = 20000.0u"K" | [2000.0, 40000.0]u"K"  | "Arrhenius temperature for lower boundary"
end
@mix struct Tup
    upperbound::typeof(1.0u"K") = 315.0u"K"   | [273.0, 325.0]u"K"     | "Upper boundary of tolerance range"
    arrupper::typeof(1.0u"K")   = 70000.0u"K" | [7000.0, 140000.0]u"K" | "Arrhenius temperature for upper boundary"
end

@Tbase mutable struct TempCorr <: AbstractTempCorr end
@Tbase @Tlow mutable struct TempCorrLower <: AbstractTempCorr end
@Tbase @Tlow @Tup mutable struct TempCorrLowerUpper <: AbstractTempCorr end


abstract type AbstractStateFeedback end

@label @range @with_kw struct Autophagy <: AbstractStateFeedback
    K_autophagy::typeof(1.0u"mol") = 0.000001u"mol" | [0.0000001u"mol", 0.00001u"mol"] | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end

abstract type AbstractAllometry end

@with_kw mutable struct SqrtAllometry <: AbstractAllometry
    scale::typeof(1.0u"m") = 0.1u"m"
end

abstract type AbstractScaling end

@label @range @with_kw struct KooijmanArea <: AbstractScaling
    M_Vref::typeof(1.0u"mol")     = 4.0u"mol"   | [0.4u"mol", 20.0u"mol"]    | "shoots scaling reference"
    M_Vscaling::typeof(1.0u"mol") = 400.0u"mol" | [40.0u"mol", 2000.0u"mol"] | "shoots scaling mass"
end

@label @range @with_kw mutable struct Params{A,S,Al,M}
    assimilation::A                        = C3Photosynthesis()      | _                             | _
    scaling::S                             = KooijmanArea()          | _                             | _
    allometry::Al                          = SqrtAllometry()         | _                             | _
    maturity::M                            = Maturity()              | _                             | _ 
    κsoma::Float64                         = 0.6                     | [0.0, 1.0]                    | "reserve flux allocated to growth"
    M_Vgerm::typeof(1.0u"mol")             = 0.5u"mol"               | _                             | "structural mass at germination"
    j_E_mai::typeof(1.0u"mol*mol^-1*d^-1") = 0.001u"mol*mol^-1*d^-1" | [0.0u, 0.01]u"mol*mol^-1*d^-1"| "spec somatic maint costs."
    j_P_mai::typeof(1.0u"mol*mol^-1*d^-1") = 0.01u"mol*mol^-1*d^-1"  | [0.0u, 0.1]u"mol*mol^-1*d^-1" | "product formation linked to maintenance"
    #
    y_P_V::typeof(1.0u"mol*mol^-1")        = 0.02u"mol*mol^-1"       | _                             | "product formation linked to growth" # where does this carbon come from?"
    y_V_E::typeof(1.0u"mol*mol^-1")        = 0.7u"mol*mol^-1"        | _                             | "from reserve to structure" # 0.3 lost as CO2?
    y_E_ET::typeof(1.0u"mol*mol^-1")       = 0.8u"mol*mol^-1"        | _                             | "translocated reserve:"
    y_EC_ECT::typeof(1.0u"mol*mol^-1")     = 1.0u"mol*mol^-1"        | _                             | "translocated C-reserve"
    y_EN_ENT::typeof(1.0u"mol*mol^-1")     = 1.0u"mol*mol^-1"        | _                             | "translocated N-reserve"
    #
    k_E::typeof(1.0u"mol*mol^-1*d^-1")     = 0.2u"mol*mol^-1*d^-1"   | [0.0, 1.0]u"mol*mol^-1*d^-1"  | "reserve turnover rate"
    k_EC::typeof(1.0u"mol*mol^-1*d^-1")    = 0.2u"mol*mol^-1*d^-1"   | [0.0, 1.0]u"mol*mol^-1*d^-1"  | "C-reserve turnover rate"
    k_EN::typeof(1.0u"mol*mol^-1*d^-1")    = 0.2u"mol*mol^-1*d^-1"   | [0.0, 1.0]u"mol*mol^-1*d^-1"  | "N-reserve turnover rate"
end

@label @range @with_kw mutable struct SharedParams{Fb,C}
    feedback::Fb                           = Autophagy()             | _                             | _
    tempcorr::C                            = TempCorrLowerUpper()    | _                             | _
    # These all really have to be the same for all organs
    # They probably have to be the same in P/V/M and E or the math is broken
    n_N_P::typeof(1.0u"mol*mol^-1")        = 0.0u"mol*mol^-1"        | _                             | "N/C in product (wood)"
    n_N_V::typeof(1.0u"mol*mol^-1")        = 0.15u"mol*mol^-1"       | _                             | "N/C in structure" # Shouldnt this be identical to the reserve?
    # If this > 0 it needs to be used in the model, but its not. n_N_C::typeof(1.0u"mol*mol^-1") = 0.0u"mol*mol^-1"        | _                             | "N/C in C-reserve"
    n_N_N::typeof(1.0u"mol*mol^-1")        = 10.0u"mol*mol^-1"       | _                             | "N/C in N-reserve"
    n_N_E::typeof(1.0u"mol*mol^-1")        = 0.2u"mol*mol^-1"        | _                             | "N/C in reserve" # TODO This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"

    w_P::typeof(1.0u"g*mol^-1")            = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot product (wood)"
    w_V::typeof(1.0u"g*mol^-1")            = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot structure"
    w_C::typeof(1.0u"g*mol^-1")            = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot C-reserve"
    w_N::typeof(1.0u"g*mol^-1")            = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot N-reserve"
    w_E::typeof(1.0u"g*mol^-1")            = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot reserve"

    y_E_CH_NO::typeof(1.0u"mol*mol^-1")    = 1.5u"mol*mol^-1"        | _                             | "from C-reserve to reserve, using nitrate" # 0.75 EC per E.
    y_E_EN::typeof(1.0u"mol*mol^-1")       = 0.5u"mol*mol^-1"        | _                             | "from N-reserve to reserve: 2 EN per E"
end

@label @range @with_kw mutable struct Maturity
    j_E_rep_mai::typeof(1.0u"mol*mol^-1*d^-1") = 0.001u"mol*mol^-1*d^-1" | [0.0, 0.01]u"mol*mol^-1*d^-1" | "shoots spec maturity maint costs "
    κrep::typeof(1.0)                          = 0.05                    | [0.0, 1.0]                    | "shoots reserve flux allocated to development/reprod."
    M_Vrep::typeof(1.0u"mol")                  = 10.0u"mol"              | _                             | "shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?
    w_M::typeof(1.0u"g*mol^-1")                = 25.0u"g*mol^-1"         | _                             | "mol-weight of shoot maturity reserve:"
    n_N_M::typeof(1.0u"mol*mol^-1")            = 10.0u"mol*mol^-1"       | _                             | "N/C in M-reserve"
end


@with_kw mutable struct CarbonVars
    J_L_F::typeof(1.0u"mol*m^-2*s^-1") = watts_to_light_mol(800.0)u"mol*m^-2*s^-1"     #| _ | "flux of useful photons"
    X_C::typeof(1.0u"mol*L^-1") = fraction_per_litre_gas_to_mols(400.0/1e6)u"mol*L^-1" #| _ | "carbon dioxide @ 400ppm"
    X_O::typeof(1.0u"mol*L^-1") = fraction_per_litre_gas_to_mols(0.21)u"mol*L^-1"      #| _ | "oxygen (21% volume in air) "
end

@with_kw mutable struct NitrogenVars
    # TODO work out the naming conventions here
    temp::typeof(25.0u"°C")       = 25.0u"°C"
    soilmoist::typeof(1.0)        = 1.0
    soilpot::typeof(1.0)          = 1.0
    soilpotshade::typeof(1.0)     = 1.0
    X_NH::typeof(1.0u"mmol*L^-1") = 5.0u"mmol*L^-1"  # | _ | "ammonia"
    X_NO::typeof(1.0u"mmol*L^-1") = 10.0u"mmol*L^-1" # | _ | "concentration of nitrate see e.g. [_@crawford1998molecular]"
    X_H::typeof(1.0u"mol*L^-1")   = 10.0u"mol*L^-1"  # | _ | _
end

@with_kw mutable struct Vars{V}
    assimilation::V                            = Photosynthesis.PhotoVars()
    scale::typeof(1.0)                         = 0.0
    rate::typeof(1.0u"mol*mol^-1*d^-1")        = 0.0u"mol*mol^-1*d^-1"
    θE::Float64                                = 0.0
    temp::typeof(1.0u"°C")                     = 25.0u"°C"
    k_E::typeof(1.0u"mol*mol^-1*d^-1")         = 0.2u"mol*mol^-1*d^-1"
    k_EC::typeof(1.0u"mol*mol^-1*d^-1")        = 0.2u"mol*mol^-1*d^-1"
    k_EN::typeof(1.0u"mol*mol^-1*d^-1")        = 0.2u"mol*mol^-1*d^-1"
    j_E_mai::typeof(1.0u"mol*mol^-1*d^-1")     = 0.001u"mol*mol^-1*d^-1"
    j_E_rep_mai::typeof(1.0u"mol*mol^-1*d^-1") = 0.001u"mol*mol^-1*d^-1"
    j_P_mai::typeof(1.0u"mol*mol^-1*d^-1")     = 0.01u"mol*mol^-1*d^-1"
    height::typeof(1.0u"m")                    = 0.0u"m"
end

mutable struct Organ{S,P,SH,V,F,F1,T,VR,FB,F1B}
    state::S
    name::Symbol
    params::P
    shared::SH
    vars::V
    J::F
    J1::F1
    time::T
    varsrecord::VR
    Jrecord::FB
    J1record::F1B
end

struct Organism{T,N,P}
    time::T
    shared::P
    nodes::N
end

@with_kw struct Scenario{E,T,N}
    environment::E = []
    time::T        = 0u"hr":1u"hr":1000u"hr"
    nodes::N       = (Organism(time=time),)
end
