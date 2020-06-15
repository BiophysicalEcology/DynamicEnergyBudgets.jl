using Test, DynamicEnergyBudgets, Unitful, OrdinaryDiffEq, DimensionalData

u0 = [1e-4, 1e-4, 1.0, 1e-4, 1e-4, 10.0]u"mol" # Initial value
du0 = fill(0.0u"mol/hr", 6); nothing

@testset "all state actually changes every timestep" begin
    plant = Plant(time = 0u"hr":1u"hr":8001u"hr")
    du = deepcopy(du0)
    u = deepcopy(u0)
    duref = deepcopy(du)
    plant(du, u, nothing, 1u"hr")    
    flatten(plant)
    @test all(duref .!= du)
end

@testset "diffeq works" begin
    du = deepcopy(du0)
    u = deepcopy(u0)
    environment = nothing

    plant = Plant(time=0.0u"hr":1.0u"hr":8001.0u"hr")
    plant(du, u, nothing, 1.0u"hr"); 

    @testset "unitfull state solver" begin
        prob = DiscreteProblem(plant, u0, (0.0u"hr", 1000.0u"hr"));
        @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true))
    end

    @testset "unitless state solver" begin
        prob = DiscreteProblem(plant, ustrip(u0), (0.0u"hr", 1000.0u"hr"));
        @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true))
    end

    # TODO test some actual results
end

nothing
