abstract type Assimilation end
abstract type CarbonAssimilation <: Assimilation end
abstract type NitrogenAssimilation <: Assimilation end
abstract type NH4_NO3_Assimilation <: NitrogenAssimilation end

@def kooijmanphoto begin
    k_C_binding::B = 1.0u"mol*s^-1" # :temp, "mols.s⁻¹, scaling rate for carbon dioxide"
    k_O_binding::B = 1.0u"mol*s^-1" # :temp, "mols.s⁻¹, scaling rate for oxygen"
    K_C::X = fraction_per_litre_gas_to_mols(40.0/1e6)u"mol*L^-1" # "half-saturation concentration of carbon dioxide"
    K_O::X = fraction_per_litre_gas_to_mols(0.0021)u"mol*L^-1" # "half-saturation concentration of oxygen"
    X_C::X = fraction_per_litre_gas_to_mols(400.0/1e6)u"mol*L^-1" # "carbon dioxide @ 400ppm"
    X_O::X = fraction_per_litre_gas_to_mols(0.21)u"mol*L^-1" # "oxygen (21% volume in air) "
end
 
@with_kw mutable struct KooijmanPhotosynthesis{B,X,J,A} <: CarbonAssimilation
    @kooijmanphoto
    J_L_F::J = 1.0u"mol*mol^-1*s^-1" # "mol/m²s, flux of useful photons"),
    J_L_K::J = 1.0u"mol*mol^-1*s^-1" # "mol/m²s, half-saturation flux of useful photons"),
    j_L_Amax::A = 20.0u"μmol*mol^-1*s^-1" # "umol.m⁻².s⁻¹, max spec uptake of useful photons"),
    j_C_Amax::A = 90.0u"μmol*mol^-1*s^-1" # "mol.m⁻².s⁻¹, max spec uptake of carbon dioxide"),
    j_O_Amax::A = 0.001u"μmol*mol^-1*s^-1" # "mol.m⁻².s⁻¹, max spec uptake of oxygen"),
end

@with_kw mutable struct KooijmanSLAPhotosynthesis{B,X,J,A,S} <: CarbonAssimilation
    @kooijmanphoto
    J_L_F::J = watts_to_light_mol(800.0)u"mol*m^-2*s^-1" # "mol/m²s, flux of useful photons"),
    J_L_K::J = watts_to_light_mol(300.0)u"mol*m^-2*s^-1" # "mol/m²s, half-saturation flux of useful photons"),
    j_L_Amax::A = 20.0u"μmol*m^-2*s^-1" # "umol.m⁻².s⁻¹, max spec uptake of useful photons"),
    j_C_Amax::A = 90.0u"μmol*m^-2*s^-1" # "mol.m⁻².s⁻¹, max spec uptake of carbon dioxide"),
    j_O_Amax::A = 0.001u"μmol*m^-2*s^-1" # "mol.m⁻².s⁻¹, max spec uptake of oxygen"),
    SLA::S = 9.10u"m^2*g^-1" # (5.0u"m^2*g^-1", 30.0u"m^2*g^-1"), "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"),
end

@with_kw mutable struct C3Photosynthesis{P} <: CarbonAssimilation
    photoparams::P = PlantPhysiology.PhotoParams()
end

@def nh4no3 begin
    j_NH_Amax::A = 50.0u"μmol*s^-1" # (0.1u"μmol*s^-1", 1000.0u"μmol*s^-1"), "max spec uptake of ammonia"), 
    j_NO_Amax::A = 50.0u"μmol*s^-1" # (0.1u"μmol*s^-1", 1000.0u"μmol*s^-1"), "max spec uptake of nitrate"), 
    K_NH::X = 10.0u"mmol*L^-1" # "half-saturation concentration of ammonia"),
    K_NO::X = 10.0u"mmol*L^-1" # "half-saturation concentration of nitrate"),
    K_H::X = 10.0u"mol*L^-1" # "half-saturation concentration of water"),
    X_NH::X = 5.0u"mmol*L^-1" # "ammonia"),
    X_NO::X = 10.0u"mmol*L^-1" # "concentration of nitrate see e.g. [_@crawford1998molecular]"),
    X_H::X = 10.0u"mol*L^-1" #
    ρNO::P = 0.7 # (0.0, 1.0), "weights preference for nitrate relative to ammonia."), # 1 or less but why?
    y_E_CH_NH::Y = 1.25u"mol*mol^-1" # "from roots C-reserve to reserve, using ammonia"),
end

@with_kw mutable struct Kooijman_NH4_NO3_Assimilation{A,X,P,Y} <: NH4_NO3_Assimilation
    @nh4no3
end
@with_kw mutable struct KooijmanSLA_NH4_NO3_Assimilation{A,X,P,Y,S} <: NH4_NO3_Assimilation
    @nh4no3
    SLA::S = 9.10u"m^2*g^-1" # (5.0u"m^2*g^-1", 30.0u"m^2*g^-1"), "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"),
end

abstract type AbstractStateFeedback end
@with_kw mutable struct Autophagy{A} <: AbstractStateFeedback
    K_autophagy::A = 0.000001u"mol" # (0.0000001u"mol", 0.00001u"mol")),
end

abstract type AbstractTempCorr end
@def tempcorrbase begin
    reference_temp::T = 310.0u"K" # (273.0u"K", 325.0u"K"), "temp for which rate pars are given"),
    arrh_temp::T = 2000.0u"K" # (200.0u"K", 4000.0u"K"), "Arrhenius temp"),
end
@def tempcorrlower begin
    lower_boundary::T = 280.0u"K" # (273.0u"K", 325.0u"K"), "lower boundary tolerance range"),
    arrh_lower::T = 20000.0u"K" # (2000.0u"K", 40000.0u"K"), "Arrhenius temp for lower boundary"),
end
@def tempcorrupper begin
    upper_boundary::T = 315.0u"K" # (273.0u"K", 325.0u"K"), "upper boundary tolerance range"),
    arrh_upper::T = 70000.0u"K" # (7000.0u"K", 140000.0u"K"), "Arrhenius temp for upper boundary"),
end
@with_kw mutable struct TempCorr{T} <: AbstractTempCorr
    @tempcorrbase
end
@with_kw mutable struct TempCorrLower{T} <: AbstractTempCorr
    @tempcorrbase
    @tempcorrlower
end
@with_kw mutable struct TempCorrLowerUpper{T} <: AbstractTempCorr
    @tempcorrbase
    @tempcorrlower
    @tempcorrupper
end

abstract type AbstractScaling end
@with_kw mutable struct KooijmanArea{T} <: AbstractScaling
    M_Vref::T = 4.0u"mol" # (0.4u"mol", 20.0u"mol"), (:exposed,)),
    M_Vscaling::T = 400.0u"mol" # (40.0u"mol", 2000.0u"mol"), (:exposed,), "shoots scaling mass"),
end

@with_kw mutable struct Structure
    n_N_V::typeof(0.15u"mol*mol^-1") = 0.15u"mol*mol^-1" # "N/C in structure. Shouldnt this be identical to the reserve?"),
    w_V::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot structure"),
    y_V_E::typeof(0.7u"mol*mol^-1") = 0.7u"mol*mol^-1" # (:exposed, :time), "from shoots reserve to structure: 0.3 lost as CO2?"),
end
@with_kw mutable struct Products
    j_P_mai::typeof(0.01u"mol*mol^-1*d^-1") = 0.01u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 0.1u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoot product formation linked to maintenance"),
    y_P_V::typeof(0.02u"mol*mol^-1") = 0.02u"mol*mol^-1" # (:exposed,), "shoot product formation linked to growth: where does this carbon come from?"),
    n_N_P::typeof(0.0u"mol*mol^-1") = 0.0u"mol*mol^-1" # "N/C in product (wood)"),
    w_P::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot product (wood)"),
end
@with_kw mutable struct Maturity
    n_N_M::typeof(10.0u"mol*mol^-1") = 10.0u"mol*mol^-1" # "N/C in M-reserve"),
    w_M::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot maturity reserve:"),
end
@with_kw mutable struct CarbonReserve
    κEC::typeof(0.2) = 0.2 # (0.0, 1.0), (:exposed,), "shoots non-processed C-reserve returned to C-reserve,"),
    k_EC::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots C-reserve turnover rate"),
    y_E_CH_NO::typeof(1.5u"mol*mol^-1") = 1.5u"mol*mol^-1" # (:exposed,), "from shoots C-reserve to reserve, using nitrate: 0.75 EC per E."),
    y_EC_ECT::typeof(1.0u"mol*mol^-1") = 1.0u"mol*mol^-1" # "from shoots C-reserve to roots C-reserve"),
    n_N_EC::typeof(0.0u"mol*mol^-1") = 0.0u"mol*mol^-1" # "N/C in C-reserve"),
    w_EC::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot C-reserve"),
end
@with_kw mutable struct NitrogenReserve
    κEN::typeof(0.5) = 0.5 # (0.0, 1.0), (:exposed, :time), "shoots non-processed N-reserve returned to N-reserve"),
    k_EN::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots N-reserve turnover rate"),
    y_E_EN::typeof(0.5u"mol*mol^-1") = 0.5u"mol*mol^-1" # (:exposed,), "from shoots N-reserve to reserve: 2 EN per E"),
    y_EN_ENT::typeof(1.0u"mol*mol^-1") = 1.0u"mol*mol^-1" # (:exposed,), "from roots N-reserve to shoots N-reserve"),
    n_N_EN::typeof(10.0u"mol*mol^-1") = 10.0u"mol*mol^-1" # "N/C in N-reserve"),
    w_EN::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot N-reserve"),
end
@with_kw mutable struct GeneralReserve
    k_E::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots reserve turnover rate"),
    y_E_ET::typeof(0.8u"mol*mol^-1") = 0.8u"mol*mol^-1" # (:exposed,), "from shoots reserve to roots reserve:"),
    n_N_E::typeof(0.2u"mol*mol^-1") = 0.2u"mol*mol^-1" # "N/C in reserve. This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"),
    w_E::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot reserve"),
end

@with_kw mutable struct Params{AS,SC,TC}
    assimilation::AS = C3Photosynthesis()
    scaling::SC = KooijmanArea()
    tempcorr::TC = TempCorrLowerUpper()
    # stateparams::StateParams(Structure(), Products(), Maturity(), GeneralReserve(), NitrogenReserve(), CarbonReserve())
    M_Vgerm::typeof(0.5u"mol") = 0.5u"mol" # "shoots structural mass at germination"),
    κsoma::typeof(0.6) = 0.6 # (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to soma"),
    κrep::typeof(0.05) = 0.05 # (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to development/reprod."),
    j_E_mai::typeof(0.001u"mol*mol^-1*d^-1") = 0.001u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 0.01u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots spec somatic maint costs."),
    j_E_rep_mai::typeof(0.001u"mol*mol^-1*d^-1") = 0.001u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 0.01u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots spec maturity maint costs "),
    M_Vrep::typeof(10.0u"mol") = 10.0u"mol" # "shoots structural mass at start reproduction"), # TODO: isn't this variable/seasonally triggered?
    # State dependent. Separate these out.
    n_N_V::typeof(0.15u"mol*mol^-1") = 0.15u"mol*mol^-1" # "N/C in structure. Shouldnt this be identical to the reserve?"),
    w_V::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot structure"),
    y_V_E::typeof(0.7u"mol*mol^-1") = 0.7u"mol*mol^-1" # (:exposed, :time), "from shoots reserve to structure: 0.3 lost as CO2?"),
    y_P_V::typeof(0.02u"mol*mol^-1") = 0.02u"mol*mol^-1" # (:exposed,), "shoot product formation linked to growth: where does this carbon come from?"),
    n_N_P::typeof(0.0u"mol*mol^-1") = 0.0u"mol*mol^-1" # "N/C in product (wood)"),
    w_P::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot product (wood)"),
    j_P_mai::typeof(0.01u"mol*mol^-1*d^-1") = 0.01u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 0.1u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoot product formation linked to maintenance"),
    n_N_M::typeof(10.0u"mol*mol^-1") = 10.0u"mol*mol^-1" # "N/C in M-reserve"),
    w_M::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot maturity reserve:"),
    κEC::typeof(0.2) = 0.2 # (0.0, 1.0), (:exposed,), "shoots non-processed C-reserve returned to C-reserve,"),
    k_EC::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots C-reserve turnover rate"),
    y_E_CH_NO::typeof(1.5u"mol*mol^-1") = 1.5u"mol*mol^-1" # (:exposed,), "from shoots C-reserve to reserve, using nitrate: 0.75 EC per E."),
    y_EC_ECT::typeof(1.0u"mol*mol^-1") = 1.0u"mol*mol^-1" # "from shoots C-reserve to roots C-reserve"),
    n_N_EC::typeof(0.0u"mol*mol^-1") = 0.0u"mol*mol^-1" # "N/C in C-reserve"),
    w_EC::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot C-reserve"),
    κEN::typeof(0.5) = 0.5 # (0.0, 1.0), (:exposed, :time), "shoots non-processed N-reserve returned to N-reserve"),
    k_EN::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots N-reserve turnover rate"),
    y_E_EN::typeof(0.5u"mol*mol^-1") = 0.5u"mol*mol^-1" # (:exposed,), "from shoots N-reserve to reserve: 2 EN per E"),
    y_EN_ENT::typeof(1.0u"mol*mol^-1") = 1.0u"mol*mol^-1" # (:exposed,), "from roots N-reserve to shoots N-reserve"),
    n_N_EN::typeof(10.0u"mol*mol^-1") = 10.0u"mol*mol^-1" # "N/C in N-reserve"),
    w_EN::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot N-reserve"),
    k_E::typeof(0.2u"mol*mol^-1*d^-1") = 0.2u"mol*mol^-1*d^-1" # (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots reserve turnover rate"),
    y_E_ET::typeof(0.8u"mol*mol^-1") = 0.8u"mol*mol^-1" # (:exposed,), "from shoots reserve to roots reserve:"),
    n_N_E::typeof(0.2u"mol*mol^-1") = 0.2u"mol*mol^-1" # "N/C in reserve. This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"),
    w_E::typeof(25.0u"g*mol^-1") = 25.0u"g*mol^-1" # "mol-weight of shoot reserve"),
end

mutable struct StateData{B,B2,A} <: DEDataVector{B}
    x::B2
    assim_state::A
    scaling::Float64
    rates::Float64
    J::Array{B,2}
    J1::Array{B,2}
    # Inner constructor to bulid J and J1 matching length(x)
    StateData{B,B2,A}(x::B2, assim_state::A, scaling, rates) where {B,B2,A} = begin
        J = zeros(B, length(TRANS), length(x))
        J1 = zeros(B, length(TRANS1), length(x))
        new(x, assim_state, scaling, rates, J, J1)
    end
end
# Outer construtor to switch B and B2 params from DEDataVector
StateData{B,A}(x::AbstractVector{B}, assim_state::A, args...) = begin
    B2 = typeof(x)
    StateData{B,B2,A}(x, assim_state, args...)
end
# Outer construtor for defaults, without J and J1
StateData(; x = StatePVE(), 
          assim_state = PhotoIntermediate(), 
          scaling = 0.0, 
          rates = 1.0) = begin
    StateData(x, assim_state, scaling, rates)
end

type OrganName{S} end

function default_params(::OrganName{:Root})
    Params(assimilation=Kooijman_NH4_NO3_Assimilation())
end
function default_params(::OrganName{:Shoot})
    Params()
end
function default_params(::OrganName{:Leaf})
    Params()
end

@with_kw struct Organ{B,A,P} <: AbstractMultiScaleArrayLeaf{B}
    values::StateData{B,A} = StateData() 
    name::Symbol = :Shoot
    params::P = default_params(OrganName{name}())
end

struct Organism{B,O<:Tuple{Vararg{<:Organ{B}}}} <: AbstractMultiScaleArray{B}
    nodes::O
    values::Vector{B}
    end_idxs::Vector{Int}
end
# Set kw defaults in outer constructor and construct() so end_idxs is correct
Organism(; nodes=(Organ(name=:Shoot,), Organ(name=:Root)), values=eltype(nodes[1])[]) = begin
    construct(Organism, nodes, values) 
end

mutable struct Scenario{B,E,TS,O<:Tuple{Vararg{<:Organism{B}}}} <: AbstractMultiScaleArrayHead{B}
    nodes::O
    values::Vector{B}
    end_idxs::Vector{Int}
    environment::E
    timestep::TS
end
# Set kw defaults in an outer constructor, so end_idxs is correct
Scenario(; nodes=(Organism(),), values=eltype(nodes[1])[], environment=[], timestep = 1.0u"hr") = begin
    construct(Scenario, nodes, values, environment, timestep) 
end
