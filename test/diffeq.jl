using Test, DynamicEnergyBudgets, Unitful, OrdinaryDiffEq, DimensionalData, Microclimate

u0 = [1e-4, 1e-4, 1.0, 1e-4, 1e-4, 10.0]u"mol" # Initial value
du0 = fill(0.0u"mol/hr", 6); nothing
env = MicroclimControl()

@testset "all state actually changes every timestep" begin
    plant = Plant(time = 0u"hr":1u"hr":8001u"hr"; 
        shared=SharedParams(
            resorption_pars=LosslessResorption(),
            tempcorr_pars=ParentTardieu(),
        ),
        environment=env
    );

    du = deepcopy(du0)
    u = deepcopy(u0)
    duref = deepcopy(du)
    plant(du, u, nothing, 1u"hr")    
    @test all(duref .!= du)
    @test plant.records[1].vars.rate > 0.0u"d^-1"
    @test plant.records[2].vars.rate > 0.0u"d^-1"
    @test plant.records[1].vars.scaling > 0.0
    @test plant.records[2].vars.scaling > 0.0

    @test plant.records[1].J[:C,:asi] > 0u"mol/hr"
    @test plant.records[2].J[:C,:asi] == 0u"mol/hr"
    @test plant.records[1].J[:N,:asi] == 0u"mol/hr"
    @test plant.records[2].J[:N,:asi] > 0u"mol/hr"
    @test plant.records[1].J[:V,:asi] == 0u"mol/hr"
    @test plant.records[2].J[:V,:asi] == 0u"mol/hr"

    @test plant.records[1].J[:N,:gro] < 0u"mol/hr"
    @test plant.records[2].J[:N,:gro] < 0u"mol/hr"
    @test plant.records[1].J[:C,:gro] < 0u"mol/hr"
    @test plant.records[2].J[:C,:gro] < 0u"mol/hr"
    @test plant.records[1].J[:V,:gro] > 0u"mol/hr"
    @test plant.records[2].J[:V,:gro] > 0u"mol/hr"

    @test plant.records[1].J[:C,:mai] < 0u"mol/hr"
    @test plant.records[1].J[:N,:mai] < 0u"mol/hr"
    @test plant.records[2].J[:C,:mai] < 0u"mol/hr"
    @test plant.records[2].J[:N,:mai] < 0u"mol/hr"
    @test plant.records[1].J[:V,:mai] == 0u"mol/hr"
    @test plant.records[1].J[:V,:mai] == 0u"mol/hr"

    @test plant.records[1].J[:C,:tra] == -plant.records[2].J[:C,:rej] > 0u"mol/hr"
    @test plant.records[2].J[:C,:tra] == -plant.records[1].J[:C,:rej] > 0u"mol/hr"
    @test plant.records[1].J[:N,:tra] == -plant.records[2].J[:N,:rej] > 0u"mol/hr"
    @test plant.records[2].J[:N,:tra] == -plant.records[1].J[:N,:rej] > 0u"mol/hr"
    @test plant.records[1].J[:V,:tra] == plant.records[1].J[:V,:rej] == 0u"mol/hr"
    @test plant.records[2].J[:V,:tra] == plant.records[2].J[:V,:rej] == 0u"mol/hr"

    @test plant.records[1].J[:N,:res] > 0u"mol/hr"
    @test plant.records[2].J[:N,:res] > 0u"mol/hr"
    @test plant.records[1].J[:C,:res] > 0u"mol/hr"
    @test plant.records[2].J[:C,:res] > 0u"mol/hr"
    @test plant.records[1].J[:V,:res] < 0u"mol/hr"
    @test plant.records[2].J[:V,:res] < 0u"mol/hr"

end

@testset "diffeq works" begin
    du = deepcopy(du0)
    u = deepcopy(u0)

    plant = Plant(time=0.0u"hr":1.0u"hr":8001.0u"hr"; 
        vars=(PlottableVars(), PlottableVars()),
        shared=SharedParams(
            resorption_pars=LosslessResorption(),
            tempcorr_pars=ParentTardieu(),
        ),
        environment=env
    )

    plant(du, u, nothing, 1.0u"hr"); 

    @testset "unitfull state solver" begin
        prob = DiscreteProblem(plant, u0, (0.0u"hr", 100.0u"hr"));
        @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true));
    end

    @testset "unitless state solver" begin
        prob = DiscreteProblem(plant, ustrip(u0), (0.0, 50.0));
        @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true))
    end

    plant.records[1].J.data
    plant.records[2].J
    plant.records[1].vars.rate
    plant.records[2].vars
    vJ = plant.organs[2].J

    vJ = view(plant.records[2].J, :, :, 100)

    # TODO test some actual results
end

nothing
