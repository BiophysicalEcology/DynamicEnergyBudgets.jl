@chaingang columns @label @range @units @default_kw

" Assimilation "
abstract type AbstractAssimilation end
" Paraent of all Carbon assimilation types"
abstract type AbstractCarbonAssimilation <: AbstractAssimilation end
" Paraent of all Nitrogen assimilation types"
abstract type AbstractNitrogenAssimilation <: AbstractAssimilation end
" Parent of Ammonia/Nitrate assimilation types"
abstract type AbstractNH4_NO3Assimilation <: AbstractNitrogenAssimilation end

@units @default_kw struct ConstantCarbonAssimilation{μMoMS} <: AbstractCarbonAssimilation
    uptake::μMoMS | 10.0 | u"μmol*mol^-1*s^-1"
end

@units @default_kw struct ConstantNitrogenAssimilation{μMoMS} <: AbstractNitrogenAssimilation
    uptake::μMoMS | 1.0 | u"mol*mol^-1*s^-1"
end

@mix @columns struct SLA{MG}
    SLA::MG | 9.10 | u"m^2*g^-1" | [5.0, 30.0] | "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"
end

" Uses FvCB photosynthesis model from Photosynthesis.jl "
@SLA struct C3Photosynthesis{P} <: AbstractCarbonAssimilation
    photoparams::P | Photosynthesis.EnergyBalance() | _ | _
end

" Parameters for simple photosynthesis module. With specific leaf area to convert area to mass "
@SLA struct KooijmanSLAPhotosynthesis{MoMoS,MoL,MoMoS,μMoMS,MoMS} <: AbstractCarbonAssimilation
    k_C_binding::MoMoS | 1.0 | u"mol*mol^-1*s^-1"                 | [0.1, 10.0] | "scaling rate for carbon dioxide"
    k_O_binding::MoMoS | 1.0 | u"mol*mol^-1*s^-1"                 | [0.1, 10.0] | "scaling rate for oxygen"
    K_C::MoL           | fraction_per_litre_gas_to_mols(40.0/1e6) | u"mol*L^-1"  | [0.1, 10.0] | "half-saturation concentration of carbon dioxide"
    K_O::MoL           | fraction_per_litre_gas_to_mols(0.0021) | u"mol*L^-1"    | [0.1, 10.0] | "half-saturation concentration of oxygen"
    J_L_K::MoMS        | watts_to_light_mol(300.0) | u"mol*m^-2*s^-1" | [0.1, 10.0]  | "half-saturation flux of useful photons"
    j_L_Amax::μMoMS    | 20.0  | u"μmol*m^-2*s^-1"                    | [0.1, 10.0]  | "max spec uptake of useful photons"
    j_C_Amax::μMoMS    | 90.0  | u"μmol*m^-2*s^-1"                    | [0.1, 10.0]  | "max spec uptake of carbon dioxide"
    j_O_Amax::μMoMS    | 0.001 | u"μmol*m^-2*s^-1"                    | [0.1, 10.0]  | "max spec uptake of oxygen"
end

" Parameters for Ammonia/Nitrate assimilation "
@columns struct Kooijman_NH4_NO3Assimilation{μMoMS,F,MoMo,MoL} <: AbstractNH4_NO3Assimilation
    j_NH_Amax::μMoMS | 50.0 | u"μmol*mol^-1*s^-1" | [0.1, 1000.0] | "max spec uptake of ammonia"
    j_NO_Amax::μMoMS | 50.0 | u"μmol*mol^-1*s^-1" | [0.1, 1000.0] | "max spec uptake of nitrate"
    ρNO::F           | 0.7  | _                   | [0.0, 1.0]    | "weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::MoMo  | 1.25 | u"mol*mol^-1"       | [0.0, 1.0]    | "from roots C-reserve to reserve using ammonia"
    K_NH::MoL        | 0.01 | u"mol*L^-1"         | [5.0, 20.0]   | "half-saturation concentration of ammonia"
    K_NO::MoL        | 0.01 | u"mol*L^-1"         | [5.0, 20.0]   | "half-saturation concentration of nitrate"
    K_H::MoL         | 10.0 | u"mol*L^-1"         | [5.0, 20.0]   | "half-saturation concentration of water"
end

" Parameters for lumped Nitrogen assimilation "
@columns struct N_Assimilation{μMoS,MoL} <: AbstractNitrogenAssimilation
    j_N_Amax::μMoS | 50.0 | u"μmol*mol^-1*s^-1" | [0.1, 1000.0] | "max spec uptake of ammonia"
    K_N::MoL       | 0.01 | u"mol*L^-1"         | [0.0, 1.0]    | "half-saturation concentration of nitrate"
    K_H::MoL       | 10.0 | u"mol*L^-1"         | [0.0, 1.0]    | "half-saturation concentration of water"
end

" Temperature correction parameters"
abstract type AbstractTempCorr{K} end

@mix @columns struct Tbase{K}
    reftemp::K    | 310.0   | u"K" | [273.0, 325.0]     | "Reference temperature for all rate parameters"
    arrtemp::K    | 2000.0  | u"K" | [200.0, 4000.0]    | "Arrhenius temperature"
end
@mix @columns struct Tlow{K}
    lowerbound::K | 280.0   | u"K" | [273.0, 325.0]     | "Lower boundary of tolerance range"
    arrlower::K   | 20000.0 | u"K" | [2000.0, 40000.0]  | "Arrhenius temperature for lower boundary"
end
@mix @columns struct Tup{K}
    upperbound::K | 315.0   | u"K" | [273.0, 325.0]     | "Upper boundary of tolerance range"
    arrupper::K   | 70000.0 | u"K" | [7000.0, 140000.0] | "Arrhenius temperature for upper boundary"
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
    K_autophagy::Mo | 0.000001 | u"mol" | [0.0000001, 0.00001] | "Half saturation metabolic rate for reincorporation of tissues. Necessary to not break the laws of thermodynamics!"
end

" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns struct SqrtAllometry{M} <: AbstractAllometry
    size::M | 0.1 | u"m" | [0.0, 1.0] | " allometry "
end

" Surface area scaling rules "
abstract type AbstractScaling end

" Surface areai scaling curve. Simulates growth and shade crowding later in life. "
@columns struct KooijmanArea{Mo} <: AbstractScaling
    M_Vref::Mo     | 4.0   | u"mol"   | [0.4, 20.0]    | "shoots scaling reference"
    M_Vscaling::Mo | 400.0 | u"mol"   | [40.0, 2000.0] | "shoots scaling mass"
end

" Model parameters that vary between organs "
@columns struct Params{A,S,Al,Ma,F,Mo,MoMoD,MoMo}
    assimilation::A | ConstantCarbonAssimilation() | _     | _          | _                                        
    scaling::S      | KooijmanArea()  | _                  | _          | _                                        
    allometry::Al   | SqrtAllometry() | _                  | _          | _                                        
    maturity::Ma    | Maturity()      | _                  | _          | _                                        
    κsoma::F        | 0.6             | _                  | [0.1,0.9]  | "reserve flux allocated to growth"       
    M_Vgerm::Mo     | 0.5             | u"mol"             | [0.0,1.0]  | "structural mass at germination"         
    y_P_V::MoMo     | 0.02            | u"mol*mol^-1"      | [0.0,1.0]  | "product formation linked to growth"     
    y_V_E::MoMo     | 0.7             | u"mol*mol^-1"      | [0.0,1.0]  | "from reserve to structure"              
    y_E_ET::MoMo    | 0.8             | u"mol*mol^-1"      | [0.0,1.0]  | "translocated reserve:"                  
    y_EC_ECT::MoMo  | 1.0             | u"mol*mol^-1"      | [0.0,1.0]  | "translocated C-reserve"                 
    y_EN_ENT::MoMo  | 1.0             | u"mol*mol^-1"      | [0.0,1.0]  | "translocated N-reserve"                 
    j_E_mai::MoMoD  | 0.001           | u"mol*mol^-1*d^-1" | [0.0,0.01] | "spec somatic maint costs."              
    j_P_mai::MoMoD  | 0.01            | u"mol*mol^-1*d^-1" | [0.0,0.1]  | "product formation linked to maintenance"
    k_E::MoMoD      | 0.2             | u"mol*mol^-1*d^-1" | [0.0,1.0]  | "reserve turnover rate"                  
    k_EC::MoMoD     | 0.2             | u"mol*mol^-1*d^-1" | [0.0,1.0]  | "C-reserve turnover rate"                
    k_EN::MoMoD     | 0.2             | u"mol*mol^-1*d^-1" | [0.0,1.0]  | "N-reserve turnover rate"                
end

" Model parameters shared between organs "
@columns struct SharedParams{Fb,C,MoMo,GMo}
    feedback::Fb    | nothing             | _ | _         | _
    tempcorr::C     | nothing             | _ | _         | _
    # n_N_P::MoMo   | 0.0u"mol*mol^-1"        | _         | "N/C in product (wood)"
    # n_N_V::MoMo   | 0.15u"mol*mol^-1"       | _         | "N/C in structure" # Shouldnt this be identical to the reserve?
    # n_N_C::MoMo   | 0.0u"mol*mol^-1"        | _         | "N/C in C-reserve"
    # n_N_N::MoMo   | 10.0u"mol*mol^-1"       | _         | "N/C in N-reserve"
    n_N_E::MoMo     | 0.2  | u"mol*mol^-1" | [0.0, 2.0]   | "N/C in reserve" # TODO This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"
    w_P::GMo        | 25.0 | u"g*mol^-1"   | [10.0, 40.0] | "mol-weight of shoot product (wood)"
    w_V::GMo        | 25.0 | u"g*mol^-1"   | [10.0, 40.0] | "mol-weight of shoot structure"
    w_C::GMo        | 25.0 | u"g*mol^-1"   | [10.0, 40.0] | "mol-weight of shoot C-reserve"
    w_N::GMo        | 25.0 | u"g*mol^-1"   | [10.0, 40.0] | "mol-weight of shoot N-reserve"
    w_E::GMo        | 25.0 | u"g*mol^-1"   | [10.0, 40.0] | "mol-weight of shoot reserve" 
    y_E_CH_NO::MoMo | 1.5  | u"mol*mol^-1" | [0.0, 2.0]   | "from C-reserve to reserve, using nitrate"
    y_E_EN::MoMo    | 1.5  | u"mol*mol^-1" | [0.0, 2.0]   | "from N-reserve to reserve"
end

" Maturity parameters. Seperated to make maturity modeling optional, reducing complexity "
@columns struct Maturity{MoMoD,F,Mo,GMo,MoMo}
    j_E_rep_mai::MoMoD | 0.001 | u"mol*mol^-1*d^-1" | [0.0, 0.01] | "shoots spec maturity maint costs "
    κrep::F            | 0.05  | _                  | [0.0, 1.0]  | "shoots reserve flux allocated to development/reprod."
    M_Vrep::Mo         | 10.0  | u"mol"             | [0.0, 1.0]  | "shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?
    w_M::GMo           | 25.0  | u"g*mol^-1"        | [0.0, 1.0]  | "mol-weight of shoot maturity reserve:"
    n_N_M::MoMo        | 10.0  | u"mol*mol^-1"      | [0.0, 1.0]  | "N/C in M-reserve"
end

" Variables for carbon assimilation "
@units @default_kw mutable struct CarbonVars{MoMS,MoL}
    J_L_F::MoMS | watts_to_light_mol(800.0)                 | u"mol*m^-2*s^-1" #| "flux of useful photons"
    X_C::MoL    | fraction_per_litre_gas_to_mols(400.0/1e6) | u"mol*L^-1"      #| "carbon dioxide @ 400ppm"
    X_O::MoL    | fraction_per_litre_gas_to_mols(0.21)      | u"mol*L^-1"      #| "oxygen (21% volume in air) "
end

" Variables for nitgroen assimilation "
@units @default_kw mutable struct NitrogenVars{F,MoL}
    # TODO work out the naming conventions here
    soilmoist::F    | 1.0  | _
    soilpot::F      | 1.0  | _
    soilpotshade::F | 1.0  | _
    X_NH::MoL       | 0.005| u"mol*L^-1"  # | "concentration of ammonia"
    X_NO::MoL       | 0.01 | u"mol*L^-1"  # | "concentration of nitrate see e.g. [_@crawford1998molecular]"
    X_H::MoL        | 10.0 | u"mol*L^-1"  # | 
end

" Model variables "
@units @default_kw mutable struct Vars{V,F,MoMoD,C,M}
    assimilation::V    | nothing | _
    scale::F           | 0.0     | _
    rate::MoMoD        | 0.0     | u"mol*mol^-1*d^-1"
    θE::F              | 0.0     | _
    temp::C            | 25.0    | u"°C"
    tempcorr::F        | 1.0     | _
    height::M          | 0.0     | u"m"
end

" Basic model components. For a plants, organs might be roots, stem and leaves "
@composite mutable struct Organ{S,P,SH,V,F,F1}
    state::S     | false
    name::Symbol | false
    params::P    | true
    shared::SH   | false
    vars::V      | false
    J::F         | false
    J1::F1       | false
end
"Construtor with keyword argument defaults"
Organ(; state = StatePVMCNE(),
        name = :Shoot,
        params = Params(),
        shared = SharedParams(),
        vars = Vars()
     ) = begin
    Organ(state, name, params, shared, vars)
end
Organ(state, name::Symbol, params, shared, vars) = begin
    one_flux = oneunit_flux(params, state)
    J = build_J(one_flux, state)
    J1 = build_J1(one_flux, state)
    Organ(state, name, params, shared, vars, J, J1)
end

Shoot() = Organ(params=Params(assimilation=ConstantCarbonAssimilation()), 
                              vars=Vars(assimilation=nothing))
Root() = Organ(params=Params(assimilation=ConstantNitrogenAssimilation()), 
                             vars=Vars(assimilation=nothing))

"Records of variables and flux for ploting and analysis"
struct Records{V,F,F1}
    vars::V
    J::F
    J1::F1
end
"Constructor for records. Arrays use the length of the current timespan"
Records(o::Organ, time) = begin
    varsrec = build_record(o.vars, time)
    Jrec = build_record(o.J, time)
    J1rec = build_record(o.J1, time)
    Records(varsrec, Jrec, J1rec)
end

"An organism, made up of organs"
@composite struct Organism{O,S,E,R}
    organs::O      | true
    shared::S      | true
    environment::E | false
    records::R     | false
end
Organism(organs::O, shared::S, environment::E, records::R) where {O,S,E,R} = begin
    organs = ([Organ(o.state, o.name, o.params, shared, o.vars) for o in organs]...)
    Organism{typeof(organs),S,E,R}(organs, shared, environment, records)
end
"Outer construtor for defaults"
Organism(; organs = (Shoot(), Root()),
           shared = SharedParams(),
           environment = nothing,                       
           time = 0u"hr":1u"hr":1000u"hr") = begin
    records = []
    for organ in organs
        push!(records, Records(organ, time))
    end
    Organism(organs, shared, environment, tuple(records...))
end

UntypedVars(;kwargs...) = default_kw(Vars{Any,Any,Any,Any,Any}; kwargs...)
UntypedOrgan(state::S, name::N, params::P, shared::Sh, vars::V) where {S,N,P,Sh,V} = begin
    # Get flux units from params, instead of explicitly.
    one_flux = oneunit_flux(params, state)
    J = build_J(one_flux, state; typ = Any)
    J1 = build_J1(one_flux, state; typ = Any)
    Organ{S,N,P,Sh,V,typeof(J),typeof(J1)}(state, name, params, shared, vars, J, J1)
end
UntypedShoot() = Organ(vars=UntypedVars(assimilation=nothing), 
                     params=Params(assimilation=ConstantCarbonAssimilation()))
UntypedRoot() = Organ(vars=UntypedVars(assimilation=nothing), 
                    params=Params(assimilation=ConstantNitrogenAssimilation()))
UntypedOrganism() = Organism(time=0:1:1000, organs=(UntypedShoot(), UntypedRoot())) 
