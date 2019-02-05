# Parameter helper functions
const gas_molpL = 22.4

watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L
fraction_per_litre_gas_to_mols(frac) = frac / 22.4

" Assimilation "
abstract type AbstractAssim end

" Parent of all Carbon assimilation types"
abstract type AbstractCAssim <: AbstractAssim end

" Parent of all Nitrogen assimilation types"
abstract type AbstractNAssim <: AbstractAssim end

" Parent of Ammonia/Nitrate assimilation types"
abstract type AbstractNH4_NO3Assim <: AbstractNAssim end

@columns struct ConstantCAssim{μMoMS} <: AbstractCAssim
    uptake::μMoMS | 0.1 | μmol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [0.0, 1.0] | _ | "constant rate of C uptake"
end

@columns struct ConstantNAssim{μMoMS} <: AbstractNAssim
    uptake::μMoMS | 0.03 | μmol*mol^-1*s^-1  | Gamma(2.0, 2.0) | [0.0, 0.5] | _ | "constant rate of N uptake" 
end

@mix @columns struct SLA{MG}
    # Field | Default | Unit     | Prior         | Limits       | Description
    SLA::MG | 9.10    | m^2*g^-1 | Gamma(10.0, 1.0) | [5.0, 30.0] | _ | "Specific leaf Area. Ferns: 17.4, Forbs: 26.2, Graminoids: 24.0, Shrubs: 9.10, Trees: 8.30"
end

" Uses FvCB photosynthesis model from Photosynthesis.jl "
@SLA struct FvCBPhotosynthesis{P} <: AbstractCAssim
    photoparams::P | Photosynthesis.FvCBEnergyBalance() | _ | _ | _ | _ | _
end

@SLA @mix @columns struct KooijmanPhoto{μMoMoS,MoL,μMoMS,MoMS}
    # Field             | Default           | Unit             | Prior           | Limits        | Log  | Description
    k_C_binding::μMoMoS | 10000.0           | μmol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [1e-5, 2000.0]   | true | "Scaling rate for carbon dioxide"
    k_O_binding::μMoMoS | 10000.0           | μmol*mol^-1*s^-1 | Gamma(2.0, 2.0) | [1e-5, 2000.0]   | true | "Scaling rate for oxygen"
    K_C::MoL            | 50*1e-6/gas_molpL | mol*L^-1         | Gamma(2.0, 2.0) | [1e-7, 100.0]   | true | "Half-saturation concentration of carbon dioxide"
    K_O::MoL            | 0.0021/gas_molpL  | mol*L^-1         | Gamma(2.0, 2.0) | [1e-7, 100.0]   | true | "Half-saturation concentration of oxygen"
    J_L_K::MoMS         | 2000.0            | mol*m^-2*s^-1    | Gamma(2.0, 2.0) | [1e-3, 100000.0] | true | "Half-saturation flux of useful photons"
    j_L_Amax::μMoMS     | 100.01            | μmol*m^-2*s^-1   | Gamma(2.0, 2.0) | [1e-4, 10000.0]  | true | "Max spec uptake of useful photons"
    j_C_Amax::μMoMS     | 20.0              | μmol*m^-2*s^-1   | Gamma(2.0, 2.0) | [5.0,  100.0]    | true | "Max spec uptake of carbon dioxide"
    j_O_Amax::μMoMS     | 0.1               | μmol*m^-2*s^-1   | Gamma(2.0, 2.0) | [0.01,  50.0]    | true | "Max spec uptake of oxygen"
end

abstract type AbstractKooijmanPhoto <: AbstractCAssim end

" Parameters for simple photosynthesis module. With specific leaf area to convert area to mass "
@KooijmanPhoto struct KooijmanSLAPhotosynthesis{} <: AbstractKooijmanPhoto end


@KooijmanPhoto struct KooijmanWaterPotentialPhotosynthesis{PL} <: AbstractKooijmanPhoto 
    potential_modifier::PL | Photosynthesis.ZhouPotentialDependence() | _ | _ | _ | _ | "modify photosynthesis with a water potential model"
end

# @columns struct WaterPotentialCutoff{KPA}
    # cutoff::KPA    | -250.0  | kPa  | Gamma(2.0, 2.0) | [0.0, -500.0] | _ | "Max spec uptake of oxygen"
# end


" Parameters for Ammonia/Nitrate assimilation "
@columns struct KooijmanNH4_NO3Assim{μMoMS,F,MoMo,MoL} <: AbstractNH4_NO3Assim
    #Field           | Default | Unit                | Prior            | Limits         | Description
    j_NH_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | _ | "Max spec uptake of ammonia"
    j_NO_Amax::μMoMS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | _ | "Max spec uptake of nitrate"
    ρNO::F           | 0.7     | _                | Beta(0.7, 1.0)   | [1e-4, 1.0]   | _ | "Weights preference for nitrate relative to ammonia." # 1 or less but why?
    y_E_CH_NH::MoMo  | 1.25    | mol*mol^-1       | Gamma(1.25, 1.0) | [1e-4, 4.0]   | _ | "From roots C-reserve to reserve using ammonia"
    K_NH::MoL        | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [1e-4, 10.0]  | _ | "Half-saturation concentration of ammonia"
    K_NO::MoL        | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [1e-4, 10.0]  | _ | "Half-saturation concentration of nitrate"
    K_H::MoL         | 10.0    | mol*L^-1         | Gamma(10.0, 1.0) | [5.0, 20.0]   | _ | "Half-saturation concentration of water"
end

" Parameters for lumped Nitrogen assimilation "
@columns struct NAssim{μMoS,MoL} <: AbstractNAssim
    # Field        | Default | Unit                | Prior            | Limits         | Description
    j_N_Amax::μMoS | 50.0    | μmol*mol^-1*s^-1 | Gamma(50.0, 1.0) | [0.1, 1000.0] | _ | "Max spec uptake of ammonia"
    K_N::MoL       | 0.01    | mol*L^-1         | Gamma(0.01, 1.0) | [1e-4, 1.0]    | _ | "Half-saturation concentration of nitrate"
    K_H::MoL       | 10.0    | mol*L^-1         | Gamma(10.0, 1.0) | [1e-2, 100.0]  | _ | "Half-saturation concentration of water"
end


" Variables for carbon assimilation "
@description @units @default_kw mutable struct CarbonVars{MoMS,MoL,KPA}
    J_L_F::MoMS              | watts_to_light_mol(800.0) | mol*m^-2*s^-1 | "flux of useful photons"
    X_C::MoL                 | (400.0/1e6)               | mol*L^-1      | "carbon dioxide concentration in air"
    X_O::MoL                 | 0.21 * gas_molpL          | mol*L^-1      | "oxygen concentration in air"
    soilwaterpotential::KPA  | 0.0                       | kPa           | "soil water potential"
end

" Variables for nitgroen assimilation "
@description @units @default_kw mutable struct NitrogenVars{F,KPA,MoL}
    # TODO work out the naming conventions here
    soilwaterpotential::KPA    | 1.0   | kPa      | "soil water potential"
    soilwaterconent::F         | 1.0   | _        | _
    X_NH::MoL                  | 0.005 | mol*L^-1 | "concentration of ammonia"
    X_NO::MoL                  | 0.01  | mol*L^-1 | "concentration of nitrate see e.g. [_@crawford1998molecular]"
    X_H::MoL                   | 10.0  | mol*L^-1 | _
end


"""
    assimilation!(f::AbstractAassim, o, u)
Runs nitrogen uptake, and combines N with translocated C.
"""
assimilation!(p::HasCNE, f::AbstractCAssim, o, u) = begin
    J = flux(o)
    c_uptake = photosynthesis(f, o, u) * u.V * shape(o)
    n_tra = J[:N,:tra]
    # Merge rejected N from root and photosynthesized C into reserves
    J[:C,:ass], J[:N,:tra], J[:E,:ass] = stoich_merge(c_uptake, n_tra, y_E_CH_NO(o), y_E_EN(o))

    # lc, ln = stoich_merge_losses(c_uptake, n_tra, J[:C,:ass], J[:N,:tra], J[:E,:ass], 
                                 # 1, 1, n_N_E(o)) 
    # J1[:C,:los] += lc
    # J1[:N,:los] += ln
end
assimilation!(::HasCN, f::AbstractCAssim, o, u) = flux(o)[:C,:ass] = begin
    println(shape(o))
    photosynthesis(f, o, u) * u.V * shape(o)
end

"""
    assimilation!(f::AbstractNH4_NO3Assim, o, u)
Runs nitrogen uptake for nitrate and ammonia, and combines N with translocated C.
Unused ammonia is discarded.
"""
assimilation!(p::HasCN, f::AbstractNH4_NO3Assim, o, u) = begin
    J = flux(o)
    J_N_ass, J_NO_ass, J_NH_ass = nitrogen_uptake(f, o, u) .* u.V * shape(o) 

    θNH = J_NH_ass/J_N_ass                          # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                   # Fraction of nitrate in arriving N-flux
    y_E_EC = θNH * f.y_E_CH_NH + θNO * shared(o).y_E_EC  # Yield coefficient from C-reserve to reserve

    c_tra = J[:C,:tra]

    # Merge rejected C from shoot and uptaken N into reserves
    J[:C,:tra], J[:N,:ass], J[:E,:ass] = stoich_merge(c_tra, J_N_ass, y_E_EC, 1/n_N_E(o))
    stoich_merge_losses(c_tra, J_N_ass, J[:C,:tra], J[:N,:ass], J[:E,:ass], 1, 1, n_N_E(o)) 

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    J[:N,:ass] = (J_NO_ass - θNO * n_N_E(o) * J[:E,:ass]) * 1/n_N_EN(o)
end

"""
    assimilation!(f::AbstractNAssim, o, u)
Runs nitrogen uptake, and combines N with translocated C.
"""
assimilation!(::HasCNE, f::AbstractNAssim, o, u) = begin
    J = flux(o)
    J_N_assim = nitrogen_uptake(f, o, u) * u.V * shape(o) 
    c_tra = J[:C,:tra]

    # Merge rejected C from shoot and uptaken N into reserves
    # treating N as N reserve now carbon has been incorporated.
    J[:C,:tra], J[:N,:ass], J[:E,:ass] = stoich_merge(c_tra, J_N_assim, y_E_CH_NO(o), y_E_EN(o))
    # lc, ln = stoich_merge_losses(c_tra, J_N_assim, J[:C,:tra], J[:N,:ass], J[:E,:ass], 1, 1, n_N_E(o)) 
    # J1[:C,:los] += lc
    # J1[:N,:los] += ln
end
assimilation!(::HasCN, f::AbstractNAssim, o, u) = flux(o)[:N,:ass] = nitrogen_uptake(f, o, u) * u.V * shape(vars(o))

"""
    photosynthesis(f::ConstantCAssim, o, u)
Returns a constant rate of carbon assimilation.
"""
photosynthesis(f::ConstantCAssim, o, u) = f.uptake

"""
    photosynthesis(f::FvCBPhotosynthesis, o, u)
Returns carbon assimilated in mols per time.
"""
photosynthesis(f::FvCBPhotosynthesis, o, u) = assimilation_vars(vars(o)).aleaf * f.SLA * w_V(o)

"""
    photosynthesis(f::KooijmanSLAPhotosynthesis, o, u)
Returns carbon assimilated in mols per time.
"""
photosynthesis(f::KooijmanSLAPhotosynthesis, o, u) = begin
    v = vars(o); va = assimilation_vars(v)
    mass_area_coef = w_V(o) * f.SLA
    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef

    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end

photosynthesis(f::KooijmanWaterPotentialPhotosynthesis, o, u) = begin
    v = vars(o); va = assimilation_vars(v)
    mass_area_coef = w_V(o) * f.SLA

    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef

    # Modify CO2 and O2 intake by water availability
    potential_mod = potential_dependence(f.potential_modifier, assimilation_vars(o).soilwaterpotential)
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef * potential_mod
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef * potential_mod


    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)

    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end


# photosynthesis(f::KooijmanWaterPotentialCutoffPhotosynthesis, o, u) = 
    # if assimilation_vars(o).soilwaterpotential > assimilation_pars(o).potential_cutoff
        # kooijman_photosynthesis(f, o, u)
    # else
        # zero(kooijman_photosynthesis(f, o, u))
    # end

"""
    nitrogen_uptake(f::ConstantNAssim, o, u)
Returns constant nitrogen assimilation.
"""
nitrogen_uptake(f::ConstantNAssim, o, u) = f.uptake

"""
    nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)
Returns total nitrogen, nitrate and ammonia assimilated in mols per time.
"""
function nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)
    v = vars(o); va = assimilation_vars(v)

    K1_NH = half_saturation(f.K_NH, f.K_H, va.X_H) # Ammonia saturation. va.X_H was multiplied by ox.scaling. But that makes no sense.
    K1_NO = half_saturation(f.K_NO, f.K_H, va.X_H) # Nitrate saturation
    J1_NH_ass = half_saturation(f.j_NH_Amax, K1_NH, va.X_NH) # Arriving ammonia mols.mol⁻¹.s⁻¹
    J_NO_ass = half_saturation(f.j_NO_Amax, K1_NO, va.X_NO) # Arriving nitrate mols.mol⁻¹.s⁻¹

    J_N_ass = J1_NH_ass + f.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

"""
    nitrogen_uptake(f::NAssim, o, u)
Returns nitrogen assimilated in mols per time.
"""
function nitrogen_uptake(f::NAssim, o, u)
    v = vars(o); va = assimilation_vars(v)
    # Ammonia proportion in soil water
    K1_N = half_saturation(f.K_N, f.K_H, va.X_H)
    # Arriving ammonia in mol mol^-1 s^-1
    half_saturation(f.j_N_Amax, K1_N, va.X_NO)
end

