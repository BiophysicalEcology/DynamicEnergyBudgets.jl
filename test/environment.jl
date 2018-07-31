using Revise, Unitful, DynamicEnergyBudgets

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

import DynamicEnergyBudgets: get_environment, define_organs

DynamicEnergyBudgets.get_environment(::Type{Val{:soiltemperature}}, env, h, i) = 15.4u"°C"
DynamicEnergyBudgets.get_environment(::Type{Val{:soilwatercontent}}, env, h, i) = 0.7
DynamicEnergyBudgets.get_environment(::Type{Val{:soilwaterpotential}}, env, h, i) = 100.0u"Pa"
DynamicEnergyBudgets.get_environment(::Type{Val{:airtemperature}}, env, h, i) = 23.7u"°C"
DynamicEnergyBudgets.get_environment(::Type{Val{:windspeed}}, env, h, i) = 3.7u"m*s^-1"
DynamicEnergyBudgets.get_environment(::Type{Val{:relhumidity}}, env, h, i) = 0.74
DynamicEnergyBudgets.get_environment(::Type{Val{:radiation}}, env, h, i) = 985.0u"W*m^-2"
DynamicEnergyBudgets.get_environment(::Type{Val{:par}}, env, h, i) = 3130.0u"mol*m^-2*s^-1"

u = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]

@testset "get environment" begin
    @test get_environment(Val{:par}, :not_nothing, 1, 1) == 3130.0u"mol*m^-2*s^-1"
end


@testset "apply environment to deb vars" begin
    o = DynamicEnergyBudgets.Plant();
    o1, o2 = define_organs(o, 1);
    v1 = o1.vars; v2 = o2.vars
    va1 = v1.assimilation; va2 = v2.assimilation

    apply(apply_environment!, (o1, o2), u, :not_nothing, 1)
    @test v1.temp == 23.7u"°C"
    @test va1.J_L_F == 3130.0u"mol*m^-2*s^-1"

    @test v2.temp == 15.4u"°C"
    @test upreferred(va2.X_H) ≈ 0.7 * 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"
end

@testset "apply environment to deb vars" begin
    o = DynamicEnergyBudgets.FvCBPlant();
    o1, o2 = define_organs(o, 1);
    v = o1.vars; va = v.assimilation;
    apply_environment!(o1, u, :no_env, 1)
    @test va.tair == 23.7u"°C" 
    @test va.tleaf > va.tair
    @test v.temp == va.tleaf
    @test va.windspeed == 3.7u"m*s^-1"
    @test va.rh == 0.74 
    @test va.rnet == 985.0u"W*m^-2"
    @test va.par == 3130.0u"mol*m^-2*s^-1"
    @test va.soilmoist == 0.7
    @test va.swp == 100.0u"Pa"
end

@testset "temperature correction" begin
    o = Organism()
    # @test tempcorr(80.0u"°C", o.shared.tempcorr) ≈ 0.0 atol=1e-10
    # @test tempcorr(-60.0u"°C", o.shared.tempcorr) ≈ 0.0 atol=1e-10
end
