# Parameter helper functions
const parconv = 4.57e-6
const gas_molpL = 22.4

@chain columns @description @limits @prior @units @default_kw

############################################################################################
# Assimilation Parameters

" Assimilation "
abstract type AbstractAssim end

" Parent of all Carbon assimilation types"
abstract type AbstractCAssim <: AbstractAssim end

" Parent of all Nitrogen assimilation types"
abstract type AbstractNAssim <: AbstractAssim end

" Parent of Ammonia/Nitrate assimilation types"
abstract type AbstractNH4_NO3Assim <: AbstractNAssim end

@columns struct ConstantCAssim{μMoMS} <: AbstractCAssim
    uptake::μMoMS | 1.0 | u"μmol*mol^-1*s^-1" | Gamma(2.0, 2.0) | [0.0, 20.0] | _
end

@columns struct ConstantNAssim{μMoMS} <: AbstractNAssim
    uptake::μMoMS | 0.1  | u"mol*mol^-1*s^-1"  | Gamma(2.0, 2.0) | [0.0, 2.0] | _
end

@mix @columns struct SLA{MG}
    # Field | Default | Unit        | Prior            | Limits       | Description
    SLA::MG | 9.10    | u"m^2*g^-1" | Gamma(10.0, 1.0) | [5.0, 30.0] | "Specific leaf Area. Ferns: 17.4, Forbs: 26.2, Graminoids: 24.0, Shrubs: 9.10, Trees: 8.30"
end

" Uses FvCB photosynthesis model from Photosynthesis.jl "
@SLA struct FvCBPhotosynthesis{P} <: AbstractCAssim
    photoparams::P | Photosynthesis.EnergyBalance() | _ | _ | _ | _
end

" Parameters for simple photosynthesis module. With specific leaf area to convert area to mass "
@SLA struct KooijmanSLAPhotosynthesis{MoMoS,MoL,MoMoS,μMoMS,MoMS} <: AbstractCAssim
    # Field            | Default           | Unit               | Prior           | Limits        | Description
    k_C_binding::MoMoS | 1.0               | u"mol*mol^-1*s^-1" | Gamma(2.0, 2.0) | [0.1, 10.0]  | "Scaling rate for carbon dioxide"
    k_O_binding::MoMoS | 1.0               | u"mol*mol^-1*s^-1" | Gamma(2.0, 2.0) | [0.1, 10.0]  | "Scaling rate for oxygen"
    K_C::MoL           | 50*1e-6/gas_molpL | u"mol*L^-1"        | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Half-saturation concentration of carbon dioxide"
    K_O::MoL           | 0.0021/gas_molpL  | u"mol*L^-1"        | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Half-saturation concentration of oxygen"
    J_L_K::MoMS        | 1000.0 * parconv  | u"mol*m^-2*s^-1"   | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Half-saturation flux of useful photons"
    j_L_Amax::μMoMS    | 20.0              | u"μmol*m^-2*s^-1"  | Gamma(2.0, 2.0) | [0.1, 100.0] | "Max spec uptake of useful photons"
    j_C_Amax::μMoMS    | 90.0              | u"μmol*m^-2*s^-1"  | Gamma(2.0, 2.0) | [0.1, 200.0] | "Max spec uptake of carbon dioxide"
    j_O_Amax::μMoMS    | 0.001             | u"μmol*m^-2*s^-1"  | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Max spec uptake of oxygen"
end

" Parameters for Ammonia/Nitrate assimilation "
@columns struct KooijmanNH4_NO3Assim{μMoMS,F,MoMo,MoL} <: AbstractNH4_NO3Assim
    #Field           | Default | Unit                | Prior            | Limits         | Description
    j_NH_Amax::μMoMS | 50.0    | u"μmol*mol^-1*s^-1" | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of ammonia"
    j_NO_Amax::μMoMS | 50.0    | u"μmol*mol^-1*s^-1" | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of nitrate"
    ρNO::F           | 0.7     | _                   | Beta(0.7, 1.0)   | [0.0, 1.0]    | "Weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::MoMo  | 1.25    | u"mol*mol^-1"       | Gamma(1.25, 1.0) | [0.0, 2.0]    | "From roots C-reserve to reserve using ammonia"
    K_NH::MoL        | 0.01    | u"mol*L^-1"         | Gamma(0.01, 1.0) | [5.0, 20.0]   | "Half-saturation concentration of ammonia"
    K_NO::MoL        | 0.01    | u"mol*L^-1"         | Gamma(0.01, 1.0) | [5.0, 20.0]   | "Half-saturation concentration of nitrate"
    K_H::MoL         | 10.0    | u"mol*L^-1"         | Gamma(10.0, 1.0) | [5.0, 20.0]   | "Half-saturation concentration of water"
end

" Parameters for lumped Nitrogen assimilation "
@columns struct NAssim{μMoS,MoL} <: AbstractNAssim
    # Field        | Default | Unit                | Prior            | Limits         | Description
    j_N_Amax::μMoS | 50.0    | u"μmol*mol^-1*s^-1" | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of ammonia"
    K_N::MoL       | 0.01    | u"mol*L^-1"         | Gamma(0.01, 1.0) | [0.0, 1.0]    | "Half-saturation concentration of nitrate"
    K_H::MoL       | 10.0    | u"mol*L^-1"         | Gamma(10.0, 1.0) | [0.0, 1.0]    | "Half-saturation concentration of water"
end


############################################################################################
# Other auxilary parameters

" Temperature correction parameters"
abstract type AbstractTempCorr{K} end

@mix @columns struct Tbase{K}
    # Field       | Default | Unit | Prior               | Limits              | Description
    reftemp::K    | 310.0   | u"K" | Gamma(310.0, 1.0)   | [273.0, 325.0]     | "Reference temperature for all rate parameters"
    arrtemp::K    | 2000.0  | u"K" | Gamma(2000.0, 1.0)  | [200.0, 4000.0]    | "Arrhenius temperature"
end
@mix @columns struct Tlow{K}
    lowerbound::K | 280.0   | u"K" | Gamma(280.0, 1.0)   | [273.0, 325.0]     | "Lower boundary of tolerance range"
    arrlower::K   | 20000.0 | u"K" | Gamma(20000.0, 1.0) | [2000.0, 40000.0]  | "Arrhenius temperature for lower boundary"
end
@mix @columns struct Tup{K}
    upperbound::K | 315.0   | u"K" | Gamma(315.0, 1.0)   | [273.0, 325.0]     | "Upper boundary of tolerance range"
    arrupper::K   | 70000.0 | u"K" | Gamma(70000, 1.0)   | [7000.0, 140000.0] | "Arrhenius temperature for upper boundary"
end

" Simple temperature correction parameters "
@Tbase struct TempCorr{K} <: AbstractTempCorr{K} end
" Temperature correction with lower boudn parameters"
@Tbase @Tlow struct TempCorrLower{K} <: AbstractTempCorr{K} end
" Temperature correction with lower and upper bound parameters"
@Tbase @Tlow @Tup struct TempCorrLowerUpper{K} <: AbstractTempCorr{K} end


" State feedback parameters. These modfy state based on state. "
abstract type AbstractStateFeedback end

" Autophagy. Parameters for self reabsorbtion when metabolic rates fall "
@columns struct Autophagy{Mo} <: AbstractStateFeedback
    # Field         | Default  | Unit   | Prior           | Limits                | Description
    K_autophagy::Mo | 0.000001 | u"mol" | Beta(2.0, 2.0)  | [0.0000001, 0.00001] | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end

" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns struct SqrtAllometry{M} <: AbstractAllometry
    # Field      | Default | Unit | Prior           | Limits      | Description
    allometry::M | 0.1     | u"m" | Gamma(2.0, 0.2) | [0.0, 1.0]  | "Allometric height/depth scaling"
end


" Surface area scaling rules "
abstract type AbstractScaling end

" Surface areai scaling curve. Simulates growth and shade crowding later in life. "
@columns struct KooijmanArea{Mo} <: AbstractScaling
    M_Vref::Mo     | 4.0     | u"mol" | Gamma(2.0, 0.2) | [0.4, 20.0]    | "Shoots scaling reference"
    M_Vscaling::Mo | 400.0   | u"mol" | Gamma(2.0, 0.2) | [40.0, 2000.0] | "Shoots scaling mass"
end


@columns @flattenable struct Translocation{D,P}
    destnames::D   | false | (:leaf,) | _ | _              | _         | "The organ/s translocated to"
    proportions::P | true  | (1.0,)   | _ | Beta(2.0, 2.0) | [0.0,1.0] | "The proportion of translocation sent in the first translocation. Only for inetermediaries. nothing = 100%"
end

############################################################################################
# Main DEB parameters

" Maturity parameters. Seperated to make maturity modeling optional, reducing complexity "
@columns struct Maturity{MoMoD,F,Mo}
    # Field            | Default         | Unit               | Prior           | Limits       | Description
    j_E_mat_mai::MoMoD | 0.001           | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0, 0.1]   | "Shoots spec maturity maint costs "
    κmat::F            | 0.05            | _                  | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Shoots reserve flux allocated to development/reprod."
    M_Vmat::Mo         | 10.0            | u"mol"             | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  w_M::GMo           | 25.0            | u"g*mol^-1"        | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Mol-weight of shoot maturity reserve:"
end

" Model parameters that vary between organs "
@columns struct Params{A,S,Al,Ma,Tr,F,Mo,MoMoD,MoMo}
    # Field            | Default         | Unit               | Prior           | Limits       | Description
    name::Symbol       | :organ          | _                  | _               | _            | _
    assimilation::A    | ConstantCAssim()| _                  | _               | _            | _
    scaling::S         | KooijmanArea()  | _                  | _               | _            | _
    allometry::Al      | SqrtAllometry() | _                  | _               | _            | _
    maturity::Ma       | Maturity()      | _                  | _               | _            | _
    translocation::Tr  | nothing         | _                  | _               | _            | _
    M_Vgerm::Mo        | 0.0             | u"mol"             | Gamma(2.0, 2.0) | [0.0,1.0]    | "Structural mass at germination"
    κsoma::F           | 0.6             | _                  | Beta(2.0, 2.0)  | [0.0,1.0]    | "Reserve flux allocated to growth"
    y_P_V::MoMo        | 0.02            | u"mol*mol^-1"      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Product formation linked to growth"
    y_V_E::MoMo        | 0.7             | u"mol*mol^-1"      | Beta(2.0, 2.0)  | [0.0,1.0]    | "From reserve to structure"
    y_E_ET::MoMo       | 0.8             | u"mol*mol^-1"      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated reserve:"
    y_EC_ECT::MoMo     | 1.0             | u"mol*mol^-1"      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated C-reserve"
    y_EN_ENT::MoMo     | 1.0             | u"mol*mol^-1"      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated N-reserve"
    j_E_mai::MoMoD     | 0.001           | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0,0.1]    | "Spec somatic maint costs."
    j_P_mai::MoMoD     | 0.01            | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0,0.1]    | "Product formation linked to maintenance"
    k_E::MoMoD         | 0.2             | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0,1.0]    | "Reserve turnover rate"
    k_EC::MoMoD        | 0.2             | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0,1.0]    | "C-reserve turnover rate"
    k_EN::MoMoD        | 0.2             | u"mol*mol^-1*d^-1" | Beta(2.0, 2.0)  | [0.0,1.0]    | "N-reserve turnover rate"
end

" Model parameters shared between organs "
@columns @flattenable struct SharedParams{Fb,C,MoMo,GMo}
    # Field         | _     | Default              | Unit          | Prior           | Limits       | Description
    feedback::Fb    | _     | nothing              | _             | _               | _            | _
    tempcorr::C     | _     | TempCorrLowerUpper() | _             | _               | _            | _
    n_N_P::MoMo     | _     | 0.0                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | _            | "N/C in product (wood)"
    n_N_V::MoMo     | _     | 0.2                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | _            | "N/C in structure" # Shouldnt this be identical to the reserve?
    n_N_M::MoMo     | _     | 0.2                  | u"mol*mol^-1" | Beta(2.0, 2.0)  | [0.0, 1.0]   | "N/C in M-reserve"
    n_N_E::MoMo     | _     | 0.2                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | [0.0, 1.0]   | "N/C in reserve" # TODO This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"
    n_N_EC::MoMo    | _     | 0.0                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | _            | "N/C in C-reserve"
    n_N_EN::MoMo    | _     | 10.0                 | u"mol*mol^-1" | Gamma(2.0, 2.0) | _            | "N/C in N-reserve"
    w_P::GMo        | false | 25.0                 | u"g*mol^-1"   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot product (wood)"
    w_V::GMo        | false | 25.0                 | u"g*mol^-1"   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot structure"
    w_C::GMo        | false | 25.0                 | u"g*mol^-1"   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot C-reserve"
    w_N::GMo        | false | 25.0                 | u"g*mol^-1"   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot N-reserve"
    w_E::GMo        | false | 25.0                 | u"g*mol^-1"   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot reserve"
    y_E_CH_NO::MoMo | _     | 1.5                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | [0.0, 2.0]   | "From C-reserve to reserve, using nitrate"
    y_E_EN::MoMo    | _     | 0.5                  | u"mol*mol^-1" | Gamma(2.0, 2.0) | [0.0, 2.0]   | "From N-reserve to reserve"
end



###########################################################################################
# Variables

" Variables for carbon assimilation "
@units @default_kw mutable struct CarbonVars{MoMS,MoL}
    J_L_F::MoMS | watts_to_light_mol(800.0) | u"mol*m^-2*s^-1" #| "flux of useful photons"
    X_C::MoL    | (400.0/1e6)               | u"mol*L^-1"      #| "carbon dioxide @ 400ppm"
    X_O::MoL    | 0.21 * gas_molpL          | u"mol*L^-1"      #| "oxygen (21% volume in air) "
end


get_environment(t::Type{Val{:par}}, env::M, interp, i) where M <: MicroclimateTable =
    lin_interp(env.metout, Val{:SOLR}, i) * 4.57u"mol*m^-2*s^-1"

" Variables for nitgroen assimilation "
@columns mutable struct NitrogenVars{F,MoL}
    # TODO work out the naming conventions here
    soilmoist::F    | 1.0   | _           |  _ | _ | _
    soilpot::F      | 1.0   | _           |  _ | _ | _
    soilpotshade::F | 1.0   | _           |  _ | _ | _
    X_NH::MoL       | 0.005 | u"mol*L^-1" |  _ | _ | "concentration of ammonia"
    X_NO::MoL       | 0.01  | u"mol*L^-1" |  _ | _ | "concentration of nitrate see e.g. [_@crawford1998molecular]"
    X_H::MoL        | 10.0  | u"mol*L^-1" |  _ | _ | _
end

" Model variables "
@units @default_kw mutable struct Vars{V,F,MoMoD,C,M}
    assimilation::V   | nothing | _
    scale::F          | [0.0]   | _
    rate::MoMoD       | [0.0]   | u"mol*mol^-1*d^-1"
    θE::F             | [0.0]   | _
    temp::C           | [25.0]  | u"°C"
    tempcorrection::F | [1.0]   | _
    height::M         | [0.0]   | u"m"
    t::Int            | 1       | _
end

assimilation(v) = v.assimilation
scale(v) = v.scale[v.t]
rate(v) = v.rate[v.t]
θE(v) = v.θE[v.t]
temp(v) = v.temp[v.t]
tempcorrection(v) = v.tempcorrection[v.t]
height(v) = v.height[v.t]
set_var!(v, fname, val) = getfield(v, fname)[v.t] = val

" Basic model components. For a plants, organs might be roots, stem and leaves "
struct Organ{P,S,V,F,F1}
    params::P
    shared::S
    vars::V
    J::F
    J1::F1
end

"Records of mutable variables and flux for ploting and analysis"
struct Records{V,F,F1}
    vars::V
    J::F
    J1::F1
end
"Constructor for records. Arrays use the length of the current timespan"
Records(params::Params, vars::Vars, time, val, typ) = begin
    vars = build_vars(vars, time)
    J = build_J(val, typ, time)
    J1 = build_J1(val, typ, time)
    Records(vars, J, J1)
end

"An organism, made up of organs"
@flattenable struct Organism{O,S,R,E}
    params::O      | true
    shared::S      | true
    records::R     | false
    environment::E | false
end

"Outer construtor for defaults"
Organism(; params = (ShootParams(), RootParams()),
           shared = SharedParams(),
           vars = (ShootVars(), RootVars()),
           records = nothing,
           environment = nothing,
           time = 0u"hr":1u"hr":8760u"hr") = begin
    if records == nothing
        recarray = []
        for i = 1:length(params)
            push!(recarray, Records(params[i], vars[i], time, 1.0u"mol/hr", typeof(1.0u"mol/hr")))
        end
        records = (recarray...,)
    end
    Organism(params, shared, records, environment)
end
