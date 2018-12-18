using OrdinaryDiffEq

global u0 = [0.0u"mol", 1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 1e-4u"mol", 0.0u"mol", 
      1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 10.0u"mol"] # Initial value
global du0 = fill(0.0u"mol/hr", 12); nothing

@testset "all state actually changes every timestep" begin
    global organism = Organism(time = 0u"hr":1u"hr":8001u"hr")
    global du = deepcopy(du0)
    global u = deepcopy(u0)
    global duref = deepcopy(du)
    organism(du, u, nothing, 1u"hr")    
    @test_broken all(duref .!= du) # P is 0.0. Prob because y_P_mai is broken
end

@testset "diffeq works" begin
    global du = deepcopy(du0)
    global u = deepcopy(u0)
    global environment = nothing

    global organism = Organism(time = 0u"hr":1u"hr":8001u"hr")
    organism(du, u, nothing, 1.0u"hr"); 
    global prob = DiscreteProblem(organism, u, (0u"hr", 8000u"hr"))
    global sol = solve(prob, FunctionMap(scale_by_time = true))
    global prob = DiscreteProblem(organism, u0, (0u"hr", 1000u"hr"));
    @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true)); nothing

    # TODO test some actual results
end

nothing
