using .Photosynthesis

export AbstractFvCBCAssim, BallBerryPotentialCAssim, BallBerryCAssim, EmaxCAssim 

export KooijmanWaterPotentialPhotosynthesis

abstract type AbstractFvCBCAssim <: AbstractCAssim end

@mix @flattenable @columns struct MixinFvCB{P,V,MG}
    vars::V        | false | Photosynthesis.EmaxVars() | _ | _ | _ | _
    photoparams::P | true  | nothing | _         | _           | _ | _
    SLA::MG        | true  | 24.0    | m^2*kg^-1 | (5.0, 30.0) | _ | "Specific leaf Area. Ferns: 17.4, Forbs: 26.2, Graminoids: 24.0, Shrubs: 9.10, Trees: 8.30"
end

"""
FCVB photosyntyhesis with Ball-Berry stomatal conductance

$(FIELDDOCTABLE)
"""
@MixinFvCB struct BallBerryCAssim{} <: AbstractFvCBCAssim end

@default BallBerryCAssim begin
    vars        | Photosynthesis.EmaxVars()
    photoparams | Photosynthesis.FvCBEnergyBalance(
                      photosynthesis=FvCBPhotosynthesis(
                          stomatal_conductance=BallBerryStomatalConductance(
                              soilmethod = PotentialSoilMethod()
                          ),
                          flux=Flux(),
                         ))

end

"""
FCVB photosyntyhesis with Ball-Berry stomatal conductance, 
and a soil water potential model.

$(FIELDDOCTABLE)
"""
@MixinFvCB struct BallBerryPotentialCAssim{} <: AbstractFvCBCAssim end

@default BallBerryPotentialCAssim begin
    vars        | Photosynthesis.EmaxVars()
    photoparams | Photosynthesis.FvCBEnergyBalance(
                      photosynthesis=FvCBPhotosynthesis(
                          stomatal_conductance=BallBerryStomatalConductance(
                              soilmethod = PotentialSoilMethod()
                          ),
                          flux=PotentialModifiedFlux(),
                      ))

end


"""
    photosynthesis(f::AbstractFvCBCAssimilation, o, u)

Returns carbon assimilated in mols per unit time.
"""
photosynthesis(f::AbstractFvCBCAssim, o, u) = f.vars.aleaf * f.SLA * w_V(o)


"""
    apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv)

Apply environment to shoot with `AbstractFvCBCAssim` assimilation model.
"""
apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv) = begin
    v = a.vars

    v.tair = airtemperature(shootenv)
    v.windspeed = windspeed(shootenv)
    v.rh = relhumidity(shootenv)
    v.rnet = radiation(shootenv)
    v.par = radiation(shootenv) * parconv
    v.vpd = vapour_pressure_deficit(v.tair, v.rh)
    # v.soilmoist = mean_soilwatercontent(rootenv)
    # swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
    swp = if typeof(rootenv) <: MicroclimControl
        soilwaterpotential(rootenv)
    else
        layermax(soilwaterpotential(rootenv.microclimate), rootenv)
    end
    v.swp = swp
    set_swp!(o, v.swp)

    # This runs energy balance which contains photosynthesis
    # Later in assimilation we can just use the result
    enbal!(v, assimilation_pars(o).photoparams)
    set_soilcorrection!(o, v.fsoil)

    update_temp!(o, v.tleaf)
end

"""
    KooijmanSLAPhotosynthesis(vars, k_C_binding, k_O_binding, K_C, K_O, J_L_K, j_L_Amax, j_C_Amax, j_O_Amax)

Koojman photosynthesis formulation modified by soil water potential.

This is untested and experimental.
"""
@MixinKooijmanPhoto struct KooijmanWaterPotentialPhotosynthesis{PL} <: AbstractKooijmanPhoto
    potential_modifier::PL | Photosynthesis.ZhouPotentialDependence() | _ | _ | _ | _ | "Modify photosynthesis with a water potential model"
end

photosynthesis(f::KooijmanWaterPotentialPhotosynthesis, o, u) = begin
    va = f.vars
    mass_area_coef = w_V(o) * f.SLA

    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef / 2

    # Modify CO2 and O2 intake by water availability to simulate stomatal closure.
    potentialcorrection = Photosynthesis.non_stomatal_potential_dependence(f.potential_modifier, va.soilwaterpotential)
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef * potentialcorrection * tempcorrection(o)
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef * potentialcorrection * tempcorrection(o)


    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)

    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end
