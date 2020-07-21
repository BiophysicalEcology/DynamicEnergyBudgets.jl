# Parameter helper functions
const gas_molpL = 22.4

watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L
fraction_per_litre_gas_to_mols(frac) = frac / 22.4


"""
    CarbonVars(J_L_F, X_C, X_O, soilwaterpotential)

Variables for carbon assimilation

$(FIELDDOCTABLE)
"""
@flattenable @udefault_kw @units @bounds @description mutable struct CarbonVars{MoMS,MoL,KPA}
    # Field                  | Flatn | Default i                 | Unit          | Bounds                             | Description
    J_L_F::MoMS              | false | watts_to_light_mol(800.0) | mol*m^-2*s^-1 | watts_to_light_mol.((0.0, 2000.0)) | "Flux of useful photons"
    X_C::MoL                 | false | 400.0*1e-6                | mol*L^-1      | (200.0, 700.0) .* 1e-6             | "Carbon dioxide concentration in air"
    X_O::MoL                 | false | 0.21 * gas_molpL          | mol*L^-1      | (0.1, 0.3) .* gas_molpL            | "Oxygen concentration in air"
    soilwaterpotential::KPA  | false | -100.0                    | kPa           | (0.0, -10000.0)                    | "Soil water potential"
end

"""
    NitrogenVars(soilwaterpotential, soilwaterconent, X_NH, X_NO, X_H)

Variables for nitrogen assimilation.

$(FIELDDOCTABLE)
"""
@flattenable @udefault_kw @units @bounds @description mutable struct NitrogenVars{F,KPA,MoL}
    soilwaterpotential::KPA  | false | -100.0 | kPa      | (0.0, -10000.0) | "Soil water potential"
    soilwaterconent::F       | false | -100.0 | _        | (0.0, -10000.0) | _
    X_NH::MoL                | false | 0.005  | mol*L^-1 | (0.0, 0.1)      | "Concentration of ammonia"
    X_NO::MoL                | false | 0.01   | mol*L^-1 | (0.0, 0.1)      | "Concentration of nitrate see e.g. (_@crawford1998molecular)"
    X_H::MoL                 | false | 10.0   | mol*L^-1 | (0.0, 20.0)     | _
end

" Abstract supertype of all Assimilation components"
abstract type AbstractAssim end

" Abstract supertype of all Carbon assimilation types"
abstract type AbstractCAssim <: AbstractAssim end

" Abstract supertype of all Nitrogen assimilation types"
abstract type AbstractNAssim <: AbstractAssim end

" Abstract supertype of Ammonia/Nitrate assimilation types"
abstract type AbstractNH4_NO3Assim <: AbstractNAssim end

"""
    ConstantCAssim(c_uptake)

C is assimilated at a constant rate, without control from the environment

$(FIELDDOCTABLE)
"""
@columns struct ConstantCAssim{μMoMS} <: AbstractCAssim
    # Field         | Def | Unit              | Bounds      | Log | Description
    c_uptake::μMoMS | 0.1 | μmol*mol^-1*s^-1  | (0.0, 10.0) | _   | "Constant rate of C uptake"
end

"""
    ConstantNAssim(n_uptake)

N is assimilated at a constant rate, without control from the environment

$(FIELDDOCTABLE)
"""
@columns struct ConstantNAssim{μMoMS} <: AbstractNAssim
    n_uptake::μMoMS | 0.1  | μmol*mol^-1*s^-1 | (0.0, 0.5)  | _   | "Constant rate of N uptake"
end

@mix @columns struct SLA{MG}
    SLA::MG         | 24.0 | m^2*kg^-1        | (5.0, 30.0) | _   | "Specific leaf Area. Ferns: 17.4, Forbs: 26.2, Graminoids: 24.0, Shrubs: 9.10, Trees: 8.30"
end

# Mixin paramters for multiple similar formulations
@mix @columns @SLA struct MixinKooijmanPhoto{μMoMoS,MoL,μMoMS,MoMS,V}
    #Field              | Default           | Unit             | Bounds           | Log  | Description
    vars::V             | CarbonVars()      | _                | _                | _    | _
    k_C_binding::μMoMoS | 10000.0           | μmol*mol^-1*s^-1 | (1e-5, 2000.0)   | true | "Scaling rate for carbon dioxide"
    k_O_binding::μMoMoS | 10000.0           | μmol*mol^-1*s^-1 | (1e-5, 2000.0)   | true | "Scaling rate for oxygen"
    K_C::MoL            | 50*1e-6/gas_molpL | mol*L^-1         | (1e-7, 100.0)    | true | "Half-saturation concentration of carbon dioxide"
    K_O::MoL            | 0.0021/gas_molpL  | mol*L^-1         | (1e-7, 100.0)    | true | "Half-saturation concentration of oxygen"
    J_L_K::MoMS         | 2000.0            | mol*m^-2*s^-1    | (1e-3, 100000.0) | true | "Half-saturation flux of useful photons"
    j_L_Amax::μMoMS     | 100.01            | μmol*m^-2*s^-1   | (1e-4, 10000.0)  | true | "Maximum specific uptake of useful photons"
    j_C_Amax::μMoMS     | 20.0              | μmol*m^-2*s^-1   | (5.0,  100.0)    | true | "Maximum specific uptake of carbon dioxide"
    j_O_Amax::μMoMS     | 0.1               | μmol*m^-2*s^-1   | (0.01,  50.0)    | true | "Maximum specific uptake of oxygen"
end

abstract type AbstractKooijmanPhoto <: AbstractCAssim end

"""
    KooijmanSLAPhotosynthesis(vars, k_C_binding, k_O_binding, K_C, K_O, J_L_K, j_L_Amax, j_C_Amax, j_O_Amax)

Parameters for simple photosynthesis module. With specific leaf area to convert area to mass

$(FIELDDOCTABLE)
"""
@MixinKooijmanPhoto struct KooijmanSLAPhotosynthesis{} <: AbstractKooijmanPhoto end

# @columns struct WaterPotentialCutoff{KPA}
    # cutoff::KPA    | -250.0  | kPa  | (0.0, -500.0) | _ | "Max spec uptake of oxygen"
# end


""" 
Parameters for Ammonia/Nitrate assimilation 

$(FIELDDOCTABLE)
"""
@columns struct KooijmanNH4_NO3Assim{μMoMS,F,MoMo,MoL,V} <: AbstractNH4_NO3Assim
    vars::V          | NitrogenVars() | _         | _             | _ | _
    j_NH_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | (0.1, 1000.0) | _ | "Max spec uptake of ammonia"
    j_NO_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | (0.1, 1000.0) | _ | "Max spec uptake of nitrate"
    ρNO::F           | 0.7     | _                | (1e-4, 1.0)   | _ | "Weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::MoMo  | 1.25    | mol*mol^-1       | (1e-4, 4.0)   | _ | "From roots C-reserve to reserve using ammonia"
    K_NH::MoL        | 0.01    | mol*L^-1         | (1e-4, 10.0)  | _ | "Half-saturation concentration of ammonia"
    K_NO::MoL        | 0.01    | mol*L^-1         | (1e-4, 10.0)  | _ | "Half-saturation concentration of nitrate"
    K_H::MoL         | 10.0    | mol*L^-1         | (5.0, 20.0)   | _ | "Half-saturation concentration of water"
end

""" 
Parameters for Nitrogen assimilation 

$(FIELDDOCTABLE)
"""
@columns struct NAssim{μMoS,MoL,V} <: AbstractNAssim
    vars::V          | NitrogenVars() | _         | _             | _ | _
    j_N_Amax::μMoS   | 50.0    | μmol*mol^-1*s^-1 | (0.1, 1000.0) | _ | "Max spec uptake of ammonia"
    K_N::MoL         | 0.01    | mol*L^-1         | (1e-4, 1.0)   | _ | "Half-saturation concentration of nitrate"
    K_H::MoL         | 10.0    | mol*L^-1         | (1e-2, 100.0) | _ | "Half-saturation concentration of water"
end



"""
    assimilation!(o::AbstractOrgan, u)

Runs assimilation methods, depending on formulation and state `u`.
"""
function assimilation! end
assimilation!(organs::Tuple, u) = map(assimilation!, organs, u)
assimilation!(o, u) = begin
    isgerminated(o, u) && assimilation!(has_reserves(o), assimilation_pars(o), o, u)
    nothing
end
assimilation!(::Nothing, x, o, u, env) = nothing
assimilation!(::Any, f::AbstractCAssim, o, u) = begin
    c = photosynthesis(f, o, u)
    flux(o)[:C,:asi] = max(c, zero(c)) * u[:V] * scaling(o)
end
assimilation!(p::HasCN, f::AbstractNH4_NO3Assim, o, u) = begin
    J = flux(o)
    J_N_ass, J_NO_ass, J_NH_ass = nitrogen_uptake(f, o, u) .* u[:V] * scaling(o)

    θNH = J_NH_ass/J_N_ass                          # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                   # Fraction of nitrate in arriving N-flux
    y_E_EC = θNH * f.y_E_CH_NH + θNO * shared(o).y_E_EC  # Yield coefficient from C-reserve to reserve

    C_tra = J[:C,:tra]

    # Merge rejected C from shoot and uptaken N into reserves
    J[:C,:tra], J[:N,:asi], J[:E,:asi] = stoich_merge(C_tra, J_N_ass, y_E_EC, 1/n_N_E(o))
    # stoich_merge_losses(c_tra, J_N_ass, J[:C,:tra], J[:N,:asi], J[:E,:asi], 1, 1, n_N_E(o))

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    J[:N,:asi] = (J_NO_ass - θNO * n_N_E(o) * J[:E,:asi]) * 1/n_N_EN(o)
end
"""
    assimilation!(f::AbstractNAssim, o, u)

Runs nitrogen uptake, and combines N with translocated C.
"""
assimilation!(::Any, f::AbstractNAssim, o, u) = begin
    n = nitrogen_uptake(f, o, u)
    flux(o)[:N,:asi] = max(n, zero(n)) * u[:V] * scaling(o)
end



"""
    photosynthesis(f::AbstractCAssim, o, u)

Return the assimilated C for `Organ` `o` with state variables `u`.
"""
function photosynthesis end

photosynthesis(f::ConstantCAssim, o, u) = f.c_uptake

photosynthesis(f::KooijmanSLAPhotosynthesis, o, u) = begin
    v = vars(o); va = f.vars
    mass_area_coef = w_V(o) * f.SLA

    # Photon flux is not temperature dependent, but O and C flux is.
    # See "Stylized facts in microalgal growth" Lorena, Marques, Kooijman, & Sousa.
    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef / 2
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef * tempcorrection(o)
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef * tempcorrection(o)

    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end




"""
    nitrogen_uptake(f::ConstantNAssim, o, u)

Returns constant nitrogen assimilation.
"""
nitrogen_uptake(f::ConstantNAssim, o, u) = f.n_uptake * tempcorrection(o)

"""
    nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)

Returns total nitrogen, nitrate and ammonia assimilated in mols per time.
"""
function nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)
    va = f.vars

    K1_NH = half_saturation(f.K_NH, f.K_H, va.X_H) # Ammonia saturation. va.X_H was multiplied by ox.scaling. But that makes no sense.
    K1_NO = half_saturation(f.K_NO, f.K_H, va.X_H) # Nitrate saturation
    J1_NH_ass = half_saturation(f.j_NH_Amax, K1_NH, va.X_NH) * tempcorrection(o) # Arriving ammonia mols.mol⁻¹.s⁻¹
    J_NO_ass = half_saturation(f.j_NO_Amax, K1_NO, va.X_NO) * tempcorrection(o) # Arriving nitrate mols.mol⁻¹.s⁻¹

    J_N_ass = J1_NH_ass + f.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

"""
    nitrogen_uptake(f::NAssim, o, u)

Returns nitrogen assimilated in mols per time.
"""
function nitrogen_uptake(f::NAssim, o, u)
    va = f.vars
    # Ammonia proportion in soil water
    K1_N = half_saturation(f.K_N, f.K_H, va.X_H)
    # Arriving ammonia in mol mol^-1 s^-1
    half_saturation(f.j_N_Amax, K1_N, va.X_NO) * tempcorrection(o)
end

