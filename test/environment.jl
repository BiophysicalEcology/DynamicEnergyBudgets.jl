using Test, DynamicEnergyBudgets, Microclimate, Unitful, Photosynthesis

import DynamicEnergyBudgets: get_environment, temp


rad = [1000.0, 900.0]u"W*m^-2"
snow = [1.0, 0.9]u"m"
airt = [25.0 24.0
        30.0 29.0]u"°C" .|> u"K"
ws = [1.0 2.0
      3.0 4.0]u"m*s^-1"
rh = [0.7 0.5
      0.8 0.6]
soilt = [20.0 19.0 18.0 17.0 16.0 15.0 14.0 13.0
         12.0 11.0 10.0 9.0  8.0  7.0  6.0  5.0]u"°C" .|> u"K"
soilwp = [-100.0 -90.0 -120.0 -80.0 -110.0 -70.0 -130.0 -60.0
          -200.0 -190.0 -220.0 -280.0 -110.0 -170.0 -230.0 -160.0]u"kPa"
soilwc = [0.3 0.3 0.2 0.3 0.2 0.3 0.2 0.3
          0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5]

environment = MicroclimPoint{0,Microclimate.LAYERINCREMENTS,Microclimate.NEXTLAYER,Microclimate.LAYERRANGE}(
    rad, snow, airt, rh, ws, soilt, soilwp, soilwc
)

u = [9.0,8.0,7.0]u"mol"


@testset "apply environment to deb vars" begin
    o = DynamicEnergyBudgets.Plant(environment=environment)
    o1, o2 = o.organs
    v1 = o1.vars; v2 = o2.vars
    p1 = o1.params; p2 = o2.params
    av1 = p1.assimilation_pars.vars; av2 = p2.assimilation_pars.vars

    av1.soilwaterpotential = 0.0u"kPa"
    av2.soilwaterpotential = 0.0u"kPa"
    av1.J_L_F == 0.0u"W*m^-2"

    t = 0u"hr"
    apply_environment!(o, o.organs, u, t)
    @test temp(v1) == 25.0u"°C"
    @test temp(v2) == 20.0u"°C"
    @test av1.J_L_F == 1000.0u"W*m^-2" * 4.57u"mol*W^-1*s^-1"
    @test av1.soilwaterpotential == -100.0u"kPa"
    # @test av2.soilwaterpotential == -100.0u"kPa"

    t = 1u"hr"
    apply_environment!(o, o.organs, u, t)
    @test temp(v1) == 30.0u"°C"
    @test temp(v2) == 12.0u"°C"
    @test av1.J_L_F == 900.0u"W*m^-2" * 4.57u"mol*W^-1*s^-1"
    @test av1.soilwaterpotential == -200.0u"kPa"
    # @test av2.soilwaterpotential == -200.0u"kPa"
end

@testset "apply environment to deb vars" begin
    o = DynamicEnergyBudgets.Plant(
        params=(Params(assimilation_pars=BallBerryPotentialCAssim()), Params()),
        environment=environment
    )

    o1, o2 = o.organs;
    v1 = o1.vars;
    p1 = o1.params;
    av = p1.assimilation_pars.vars

    av.tair = 0.0u"°C"
    av.windspeed = 0.0u"m*s^-1"
    av.rh = 0.0
    av.rnet = 0.0u"W*m^-2"
    av.par = 0.0u"W*m^-2" * 0.0u"mol*W^-1*s^-1"
    av.soilmoist = 0.0
    av.swp = 0.0u"kPa"

    t = 0u"hr"
    apply_environment!(o, o.organs, u, t)
    @test av.tair == 25.0u"°C"
    @test av.tleaf > av.tair
    @test temp(v1) == av.tleaf
    @test av.windspeed == 1.0u"m*s^-1"
    @test av.rh == 0.7
    @test av.rnet == 1000.0u"W*m^-2"
    @test av.par == 1000.0u"W*m^-2" * 4.57u"mol*W^-1*s^-1"
    # @test av.soilmoist == 0.2 not used currently
    @test av.swp == -100.0u"kPa"

    t = 1u"hr"
    apply_environment!(o, o.organs, u, t)
    @test av.tair == 30.0u"°C"
    @test av.tleaf > av.tair
    @test temp(v1) == av.tleaf
    @test av.windspeed == 3.0u"m*s^-1"
    @test av.rh == 0.8
    @test av.rnet == 900.0u"W*m^-2"
    @test av.par == 900.0u"W*m^-2" * 4.57u"mol*W^-1*s^-1"
    # @test av.soilmoist == 0.4 not used currently
    @test av.swp == -200.0u"kPa"
end

nothing
