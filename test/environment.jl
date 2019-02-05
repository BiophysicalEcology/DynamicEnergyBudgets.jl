import DynamicEnergyBudgets: get_environment, temp

soiltemperature = [15.4u"°C"]
soilwatercontent = [0.7]
soilwaterpotential = [100.0u"Pa"]
airtemperature = [23.7u"°C"]
windspeed = [3.7u"m*s^-1"]
relhumidity = [0.74]
radiation = [985.0u"W*m^-2"]
par = [3130.0u"mol*m^-2*s^-1"]

environment = Env(soiltemperature, soilwatercontent, soilwaterpotential, 
                  airtemperature, windspeed, relhumidity, radiation, par)

u = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]


@testset "apply environment to deb vars" begin
    global o = DynamicEnergyBudgets.Plant(environment=environment)
    global o1, o2 = define_organs(o, 1)
    global v1 = o1.vars; v2 = o2.vars
    global va1 = v1.assimilation 
    global va2 = v2.assimilation

    apply(apply_environment!, (o1, o2), u, :not_nothing, 1)
    @test temp(v1) == 23.7u"°C"
    @test va1.J_L_F == 3130.0u"mol*m^-2*s^-1"

    @test temp(v2) == 15.4u"°C"
    @test upreferred(va2.X_H) ≈ 0.7 * 1.0u"m^3*m^-3" * 1u"kg*L^-1" / 18.0u"g*mol^-1"
end

@testset "apply environment to deb vars" begin
    global o = DynamicEnergyBudgets.FvCBPlant(environment=environment)
    global o1, o2 = define_organs(o, 1)
    global v = o1.vars; 
    global va = v.assimilation;
    apply_environment!(o1, u, :no_env, 1)
    @test va.tair == 23.7u"°C" 
    @test va.tleaf > va.tair
    @test temp(v) == va.tleaf
    @test va.windspeed == 3.7u"m*s^-1"
    @test va.rh == 0.74 
    @test va.rnet == 985.0u"W*m^-2"
    @test va.par == 3130.0u"mol*m^-2*s^-1"
    @test va.soilmoist == 0.7
    @test va.swp == 100.0u"Pa"
end

nothing
