using OrdinaryDiffEq

@testset "diffeq works" begin
    environment = nothing
    u0 = [0.0u"mol", 1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 1e-4u"mol", 0.0u"mol", 
          1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 10.0u"mol"] # Initial value
    u = u0
    du = fill(0.0u"mol/hr", 12); nothing

    organism = Organism(time = 0u"hr":1u"hr":8001u"hr")
    organism(du, u0, nothing, 1.0u"hr"); 
    prob = DiscreteProblem(organism, u, (0u"hr", 8000u"hr"))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    prob = DiscreteProblem(organism, u0, (0u"hr", 1000u"hr"));
    @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true)); nothing

    # TODO test some actual results
end
