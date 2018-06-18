using Revise
using Unitful
using DynamicEnergyBudgets

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

import DynamicEnergyBudgets.get_environment

DynamicEnergyBudgets.get_environment(::Type{Val{:soiltemperature}}, env, h, i) = 15.4u"°C"
DynamicEnergyBudgets.get_environment(::Type{Val{:soilwatercontent}}, env, h, i) = 0.7
DynamicEnergyBudgets.get_environment(::Type{Val{:soilwaterpotential}}, env, h, i) = 100.0u"Pa"
DynamicEnergyBudgets.get_environment(::Type{Val{:airtemperature}}, env, h, i) = 23.7u"°C"
DynamicEnergyBudgets.get_environment(::Type{Val{:windspeed}}, env, h, i) = 3.7u"m*s^-1"
DynamicEnergyBudgets.get_environment(::Type{Val{:relhumidity}}, env, h, i) = 0.74
DynamicEnergyBudgets.get_environment(::Type{Val{:radiation}}, env, h, i) = 985.0u"W*m^-2"
DynamicEnergyBudgets.get_environment(::Type{Val{:par}}, env, h, i) = 3130.0u"mol*m^-2*s^-1"

@testset "apply environment to deb vars" begin
    o1 = Organ(params=Params(assimilation=KooijmanSLAPhotosynthesis()), vars=Vars(assimilation=CarbonVars()));
    o2 = Organ(params=Params(assimilation=N_Assimilation()), vars=Vars(assimilation=NitrogenVars()));
    va1 = o1.vars.assimilation;
    va2 = o2.vars.assimilation;

    apply(apply_environment!, (o1, o2), :no_env, 1)
    @test va1.tair == 23.7u"°C"
    @test va1.J_L_F == 3130.0u"mol*m^-2*s^-1"

    @test va2.temp == 15.4u"°C"
    @test upreferred(va2.X_H) ≈ 0.7 * 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"
end

@testset "apply environment to deb vars" begin
    o = Organ(params=Params(assimilation=C3Photosynthesis()), vars=Vars(assimilation=Photosynthesis.PhotoVars()));
    va = o.vars.assimilation;
    apply_environment!(o, :no_env, 1)
    @test va.tair == 23.7u"°C"
    @test va.windspeed == 3.7u"m*s^-1"
    @test va.rh == 0.74 
    @test va.rnet == 985.0u"W*m^-2"
    @test va.par == 3130.0u"mol*m^-2*s^-1"
    @test va.soilmoist == 0.7
    @test va.swp == 100.0u"Pa"
end

@testset "temperature correction" begin
    o = Organism()
    @test tempcorr(80.0u"°C", o.shared.tempcorr) ≈ 0.0 atol=1e-10
    @test tempcorr(-60.0u"°C", o.shared.tempcorr) ≈ 0.0 atol=1e-10
end
