# Parameter helper functions
const parconv = 4.57e-6
const gas_molpL = 22.4

watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L
fraction_per_litre_gas_to_mols(frac) = frac / 22.4

@chain columns @description @limits @prior @units @default_kw

############################################################################################
# SUs
abstract type AbstractSU end

struct ParallelComplementarySU <: AbstractSU end

struct MinimumRuleSU <: AbstractSU end

@columns struct KfamilySU{K} <: AbstractSU 
    k::K | 1.0 | _ | Gamma(2.0, 2.0) | [0.0, 10.0]  | "Synthesizing unit efficiency"
end



" Surface area scaling rules "
abstract type AbstractShape end

struct Isomorph <: AbstractShape end 
@columns struct V0morph{Mo} <: AbstractShape
    Vd::Mo  | 4.0 | mol | Gamma(2.0, 0.2) | [0.0, 1000.0] | "reference"
end
@columns struct V1morph{Mo} <: AbstractShape
    Vd::Mo  | 4.0 | mol | Gamma(2.0, 0.2) | [0.0, 1000.0] | "reference"
end
@columns struct V1V0morph{Mo,B} <: AbstractShape
    Vd::Mo    | 4.0 | mol | Gamma(2.0, 0.2) | [0.0, 1000.0] | "reference"
    Vmax::Mo  | 4.0 | mol | Gamma(2.0, 0.2) | [0.0, 1000.0] | "reference"
    β::B      | 4.0 | mol | Gamma(2.0, 0.2) | [0.0, 10.0] | "reference"
end
@columns struct Plantmorph{Mo} <: AbstractShape
    M_Vref::Mo     | 4.0     | mol | Gamma(2.0, 0.2) | [0.4, 20.0]    | "Scaling reference"
    M_Vscaling::Mo | 400.0   | mol | Gamma(2.0, 0.2) | [40.0, 2000.0] | "Scaling mass"
end

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
    uptake::μMoMS | 1.0 | μmol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [0.0, 2.0] | _
end

@columns struct ConstantNAssim{μMoMS} <: AbstractNAssim
    uptake::μMoMS | 1.0  | μmol*mol^-1*s^-1  | Gamma(2.0, 2.0) | [0.0, 2.0] | _
end

@mix @columns struct SLA{MG}
    # Field | Default | Unit        | Prior            | Limits       | Description
    SLA::MG | 9.10    | m^2*g^-1 | Gamma(10.0, 1.0) | [5.0, 30.0] | "Specific leaf Area. Ferns: 17.4, Forbs: 26.2, Graminoids: 24.0, Shrubs: 9.10, Trees: 8.30"
end

" Uses FvCB photosynthesis model from Photosynthesis.jl "
@SLA struct FvCBPhotosynthesis{P} <: AbstractCAssim
    photoparams::P | Photosynthesis.EnergyBalance() | _ | _ | _ | _
end

" Parameters for simple photosynthesis module. With specific leaf area to convert area to mass "
@SLA struct KooijmanSLAPhotosynthesis{MoMoS,MoL,MoMoS,μMoMS,MoMS} <: AbstractCAssim
    # Field            | Default           | Unit               | Prior           | Limits        | Description
    k_C_binding::MoMoS | 1.0               | mol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [0.1, 10.0]  | "Scaling rate for carbon dioxide"
    k_O_binding::MoMoS | 1.0               | mol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [0.1, 10.0]  | "Scaling rate for oxygen"
    K_C::MoL           | 50*1e-6/gas_molpL | mol*L^-1        | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Half-saturation concentration of carbon dioxide"
    K_O::MoL           | 0.0021/gas_molpL  | mol*L^-1        | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Half-saturation concentration of oxygen"
    J_L_K::MoMS        | 1000.0 * parconv  | mol*m^-2*s^-1   | Gamma(2.0, 2.0) | [0.0, 1000.0]   | "Half-saturation flux of useful photons"
    j_L_Amax::μMoMS    | 20.0              | μmol*m^-2*s^-1  | Gamma(2.0, 2.0) | [0.1, 100.0] | "Max spec uptake of useful photons"
    j_C_Amax::μMoMS    | 90.0              | μmol*m^-2*s^-1  | Gamma(2.0, 2.0) | [0.1, 200.0] | "Max spec uptake of carbon dioxide"
    j_O_Amax::μMoMS    | 0.001             | μmol*m^-2*s^-1  | Gamma(2.0, 2.0) | [0.0, 1.0]   | "Max spec uptake of oxygen"
end

" Parameters for Ammonia/Nitrate assimilation "
@columns struct KooijmanNH4_NO3Assim{μMoMS,F,MoMo,MoL} <: AbstractNH4_NO3Assim
    #Field           | Default | Unit                | Prior            | Limits         | Description
    j_NH_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of ammonia"
    j_NO_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of nitrate"
    ρNO::F           | 0.7     | _                | Beta(0.7, 1.0)   | [0.0, 1.0]    | "Weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::MoMo  | 1.25    | mol*mol^-1       | Gamma(1.25, 1.0) | [0.0, 4.0]    | "From roots C-reserve to reserve using ammonia"
    K_NH::MoL        | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [0.0, 10.0]   | "Half-saturation concentration of ammonia"
    K_NO::MoL        | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [0.0, 10.0]   | "Half-saturation concentration of nitrate"
    K_H::MoL         | 10.0    | mol*L^-1         | Gamma(10.0, 1.0) | [5.0, 20.0]   | "Half-saturation concentration of water"
end

" Parameters for lumped Nitrogen assimilation "
@columns struct NAssim{μMoS,MoL} <: AbstractNAssim
    # Field        | Default | Unit                | Prior            | Limits         | Description
    j_N_Amax::μMoS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | "Max spec uptake of ammonia"
    K_N::MoL       | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [0.0, 1.0]    | "Half-saturation concentration of nitrate"
    K_H::MoL       | 10.0    | mol*L^-1         | Gamma(10.0, 1.0) | [0.0, 100.0]  | "Half-saturation concentration of water"
end


############################################################################################
# Other auxilary parameters

" Temperature correction parameters"
abstract type AbstractTempCorr{T} end

# Temperature mixins
@mix @columns struct Tbase{T}
    # Field       | Default | Unit | Prior               | Limits              | Description
    reftemp::T    | 310.0   | K | Gamma(310.0, 1.0)   | [273.0, 325.0]     | "Reference temperature for all rate parameters"
    arrtemp::T    | 2000.0  | K | Gamma(2000.0, 1.0)  | [200.0, 4000.0]    | "Arrhenius temperature"
end
@mix @columns struct Tlow{T}
    lowerbound::T | 280.0   | K | Gamma(280.0, 1.0)   | [273.0, 325.0]     | "Lower boundary of tolerance range"
    arrlower::T   | 20000.0 | K | Gamma(20000.0, 1.0) | [2000.0, 40000.0]  | "Arrhenius temperature for lower boundary"
end
@mix @columns struct Tup{T}
    upperbound::T | 315.0   | K | Gamma(315.0, 1.0)   | [273.0, 325.0]     | "Upper boundary of tolerance range"
    arrupper::T   | 70000.0 | K | Gamma(70000, 1.0)   | [7000.0, 140000.0] | "Arrhenius temperature for upper boundary"
end

# Temperature types
" Simple temperature correction parameters "
@Tbase struct TempCorr{T} <: AbstractTempCorr{T} end
" Temperature correction with lower boudn parameters"
@Tbase @Tlow struct TempCorrLower{T} <: AbstractTempCorr{T} end
" Temperature correction with lower and upper bound parameters"
@Tbase @Tlow @Tup struct TempCorrLowerUpper{T} <: AbstractTempCorr{T} end


" State feedback parameters. These modfy state based on state. "
abstract type AbstractStateFeedback end

" Autophagy. Parameters for self reabsorbtion when metabolic rates fall "
@columns struct Autophagy{Mo} <: AbstractStateFeedback
    # Field         | Default  | Unit | Prior           | Limits               | Description
    K_autophagy::Mo | 0.000001 | mol  | Beta(2.0, 2.0)  | [0.0000001, 0.00001] | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end

" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns struct SqrtAllometry{M} <: AbstractAllometry
    # Field         | Default  | Unit | Prior           | Limits         | Description
    allometry::M    | 0.1      | m    | Gamma(2.0, 0.2) | [0.0, 1.0]   | "Allometric height/depth scaling"
end

abstract type AbstractTranslocation end
abstract type AbstractDissipativeTranslocation <: AbstractTranslocation end
abstract type AbstractLosslessTranslocation <: AbstractTranslocation end

# Dissipation mixins
@mix @columns @flattenable struct Trans{F}
    κtra::F        | 0.6     | _   | Beta(2.0, 2.0)  | [0.0,1.0]      | "Reserve flux allocated to translocation"
end

@mix @columns @flattenable struct MultiTrans{D,P}
    destnames::D   | false | (:leaf,) | _ | _              | _         | "The organ/s translocated to"
    proportions::P | true  | (1.0,)   | _ | Beta(2.0, 2.0) | [0.0,1.0] | "The proportion of translocation sent in the first translocation. Only for inetermediaries. nothing = 100%"
end

@mix @columns @flattenable struct DissTrans{MoMo}
    y_E_ET::MoMo   | true | 0.8       | mol*mol^-1      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated reserve:"
end

# Dissipation types
@Trans @MultiTrans @DissTrans struct DissipativeMultipleTranslocation{} <: AbstractDissipativeTranslocation end
@Trans @DissTrans struct DissipativeTranslocation{} <: AbstractDissipativeTranslocation  end
@Trans @MultiTrans struct LosslessMultipleTranslocation{} <: AbstractLosslessTranslocation end
@Trans struct LosslessTranslocation{} <: AbstractLosslessTranslocation end


############################################################################################
# Main DEB parameters

abstract type AbstractMaturity end

" Maturity parameters. Seperated to make maturity modeling optional, reducing complexity "
@columns struct Maturity{MoMo,MoMoD,F,Mo} <: AbstractMaturity
    # Field            | Default         | Unit            | Prior           | Limits       | Description
    n_N_M::MoMo        | 0.1             | mol*mol^-1      | Gamma(2.0, 2.0) | [0.0, 1.0]   | "N/C use for maturity"
    j_E_mat_mai::MoMoD | 0.001           | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0, 0.1]   | "Shoots spec maturity maint costs "
    κmat::F            | 0.05            | _               | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Shoots reserve flux allocated to development/reprod."
    M_Vmat::Mo         | 10.0            | mol             | Beta(2.0, 2.0)  | [0.0, 20.0]  | "Shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  w_M::GMo           | 25.0            | g*mol^-1        | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Mol-weight of shoot maturity reserve:"
end

abstract type AbstractProduction end

@columns struct Production{MoMo,MoMoD,GMo} <: AbstractProduction
    y_P_V::MoMo        | 0.02            | mol*mol^-1      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Product formation linked to growth"
    j_P_mai::MoMoD     | 0.001           | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,0.1]    | "Product formation linked to maintenance"
    n_N_P::MoMo        | 0.1             | mol*mol^-1      | Gamma(2.0, 2.0) | [0.0, 1.0]   | "N/C in product (wood)"
    w_P::GMo           | 25.0            | g*mol^-1        | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot product (wood)"
end

abstract type AbstractRejection end

@columns struct DissipativeRejection{MoMo} <: AbstractRejection
    y_EC_ECT::MoMo       | 1.0             | mol*mol^-1      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated C-reserve"
    y_EN_ENT::MoMo       | 1.0             | mol*mol^-1      | Beta(2.0, 2.0)  | [0.0,1.0]    | "Translocated N-reserve"
end
struct LosslessRejection <: AbstractRejection end

abstract type AbstractRate end
struct SimpleRate <: AbstractRate end
struct FZeroRate <: AbstractRate end

abstract type AbstractGermination end

@columns struct Germination{Mo} <: AbstractGermination
    M_Vgerm::Mo | 0.0 | mol             | Gamma(2.0, 2.0) | [0.0,1.0]    | "Structural mass at germination"
end

abstract type AbstractTurnover end

@mix @columns struct TurnoverCN{MoMoD}
    k_EC::MoMoD | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]    | "C-reserve turnover rate"
    k_EN::MoMoD | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]    | "N-reserve turnover rate"
end

@mix @columns struct TurnoverE{MoMoD}
    k_E::MoMoD  | 0.2 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]    | "C-reserve turnover rate"
end

@TurnoverCN struct TurnoverCN{} <: AbstractTurnover end
@TurnoverE @TurnoverCN struct TurnoverCNE{} <: AbstractTurnover end

# κEC::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]    | "Non-processed C-reserve returned to C-reserve"
# κEN::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]    | "Non-processed N-reserve returned to N-reserve"

abstract type AbstractMaintenance end

@columns struct Maintenance{MoMoD} <: AbstractMaintenance
    j_E_mai::MoMoD       | 0.01            | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,0.1]    | "Spec somatic maint costs."
end

abstract type AbstractParams end

" Model parameters that vary between organs "
@default_kw struct Params{A,S,Al,Ma,Tr,Re,Ge,Mt,Pr,Tu} <: AbstractParams
    # Field              | Default
    name::Symbol         | :organ
    rate_formula         | FZeroRate()
    assimilation_pars::A | ConstantCAssim()
    shape_pars::S        | Plantmorph()
    allometry_pars::Al   | nothing
    maturity_pars::Ma    | nothing
    trans_pars::Tr       | nothing
    rejection_pars::Re   | LosslessRejection()
    germination_pars::Ge | Germination()
    maintenance_pars::Mt | Maintenance()
    production_pars::Pr  | Production()
    turnover_pars::Tu    | TurnoverCN()
end

abstract type AbstractCore end

@columns struct Core{MoMo,GMo} <: AbstractCore
    y_V_E::MoMo  | 0.7   | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0,1.0]    | "From reserve to structure"
    y_E_EC::MoMo | 1.0   | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 4.0]   | "From C-reserve to reserve, using nitrate"
    y_E_EN::MoMo | 1.0   | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 4.0]   | "From N-reserve to reserve"
    n_N_V::MoMo  | 0.15  | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 1.0]   | "N/C in structure"
    n_N_E::MoMo  | 0.2   | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 1.0]   | "N/C in reserve"
    w_V::GMo     | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot structure"
    # w_C::GMo   | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot C-reserve"
    # w_N::GMo   | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot N-reserve"
    # w_E::GMo   | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [10.0, 40.0] | "Mol-weight of shoot reserve"
end


" Model parameters shared between organs "
@default_kw struct SharedParams{SU,Co,FB,TC}
    su_pars::SU       | ParallelComplementarySU()
    core_pars::Co     | Core()
    feedback_pars::FB | nothing
    tempcorr_pars::TC | nothing
end



###########################################################################################
# Variables

" Variables for carbon assimilation "
@description @units @default_kw mutable struct CarbonVars{MoMS,MoL}
    J_L_F::MoMS | watts_to_light_mol(800.0) | mol*m^-2*s^-1 | "flux of useful photons"
    X_C::MoL    | (400.0/1e6)               | mol*L^-1      | "carbon dioxide @ 400ppm"
    X_O::MoL    | 0.21 * gas_molpL          | mol*L^-1      | "oxygen (21% volume in air) "
end

" Variables for nitgroen assimilation "
@description @units @default_kw mutable struct NitrogenVars{F,MoL}
    # TODO work out the naming conventions here
    soilmoist::F    | 1.0   | _        | _
    soilpot::F      | 1.0   | _        | _
    soilpotshade::F | 1.0   | _        | _
    X_NH::MoL       | 0.005 | mol*L^-1 | "concentration of ammonia"
    X_NO::MoL       | 0.01  | mol*L^-1 | "concentration of nitrate see e.g. [_@crawford1998molecular]"
    X_H::MoL        | 10.0  | mol*L^-1 | _
end

" Model variables "
@units @default_kw mutable struct Vars{V,F,MoMoD,C,M}
    assimilation_vars::V | nothing | _
    shape::F             | [0.0]   | _
    rate::MoMoD          | [0.0]   | mol*mol^-1*d^-1
    θE::F                | [0.0]   | _
    temp::C              | [25.0]  | °C
    tempcorrection::F    | [1.0]   | _
    height::M            | [0.0]   | m
    t::Int               | 1       | _
end

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
Records(params, vars, time, val, typ) = begin
    vars = build_vars(vars, time)
    J = build_J(val, typ, time)
    J1 = build_J1(val, typ, time)
    Records(vars, J, J1)
end

"An organism, made up of organs"
@flattenable struct Organism{P,S,R,O,E,D}
    params::P      | true
    shared::S      | true
    records::R     | false
    organs::O      | false
    environment::E | false
    dead::D        | false
end

"Outer construtor for defaults"
Organism(; params = (ShootParamsCN(), RootParamsCN()),
           shared = SharedParams(),
           vars = (ShootVars(), RootVars()),
           records = nothing,
           environment = nothing,
           time = 0hr:1hr:8760hr) = begin
    if records == nothing
        records = []
        for i = 1:length(params)
            push!(records, Records(params[i], vars[i], time, 1.0mol/hr, typeof(1.0mol/hr)))
        end
        records = (records...,)
    end
    organs = define_organs(params, shared, records, 1)
    dead = Array{Bool,0}(undef)
    dead[] = false
    Organism(params, shared, records, organs, environment, dead)
end


# Traits

struct HasCN end
struct HasCNE end

has_reserves(o) = 
  if typeof(turnover_pars(o)) <: TurnoverCN 
      HasCN()
  elseif typeof(turnover_pars(o)) <: TurnoverCNE
      HasCNE()
  else
      nothing
  end

