using OrdinaryDiffEq

@testset "diffeq works" begin
    # TODO test some actual results
    environment = nothing
    u0 = [0.0u"mol", 1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 1e-4u"mol", 0.0u"mol", 
          1e-4u"mol", 0.0u"mol", 1e-4u"mol", 1e-4u"mol", 10.0u"mol"] # Initial value
    u = u0
    du = fill(0.0u"mol/hr", 12); nothing

    organism = Organism(nodes=(Organ(params=Params(assimilation=KooijmanSLAPhotosynthesis())), 
                        Organ(params=Params(assimilation=N_Assimilation()))),
                       );
    runmodel!(du, u0, organism, 1.0u"hr"); 
    prob = DiscreteProblem(runmodel!, u0, (0u"hr", 1000u"hr"), organism);
    @test_nowarn sol = solve(prob, FunctionMap(scale_by_time = true)); nothing
end
