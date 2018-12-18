function nfactory()
    o1 = construct_organ(params=Params(y_EC_ECT=0.8, y_EN_ENT=0.8),
               vars=Vars(assimilation=NitrogenVars()));;
    p1 = o1.params
    o2 = construct_organ(params=Params(y_EC_ECT=0.8, y_EN_ENT=0.8))
    p2 = o2.params;
    u1 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u1 = LVector{eltype(u1),typeof(u1),DynamicEnergyBudgets.STATE}(u1)
    u2 = LVector{eltype(u2),typeof(u2),DynamicEnergyBudgets.STATE}(u2)
    u2 .*= 2.7
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    du1 = LVector{eltype(du1),typeof(du1),DynamicEnergyBudgets.STATE}(du1)
    du2 = LVector{eltype(du2),typeof(du2),DynamicEnergyBudgets.STATE}(du2)
    set_var!(o1.vars, :scale, scaling(o1.params.scaling, u1[:V]))
    set_var!(o2.vars, :scale, scaling(o2.params.scaling, u2[:V]))
    f = NAssim();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

function cfactory()
    o1 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8),
               vars=Vars(assimilation=CarbonVars()));;
    p1 = o1.params;
    o2 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8));
    p2 = o2.params;
    u1 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u1 = LVector{eltype(u1),typeof(u1),DynamicEnergyBudgets.STATE}(u1)
    u2 = LVector{eltype(u2),typeof(u2),DynamicEnergyBudgets.STATE}(u2)
    u2 .*= 1.9
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    du1 = LVector{eltype(du1),typeof(du1),DynamicEnergyBudgets.STATE}(du1)
    du2 = LVector{eltype(du2),typeof(du2),DynamicEnergyBudgets.STATE}(du2)
    set_var!(o1.vars, :scale, scaling(o1.params.scaling, u1[:V]))
    set_var!(o2.vars, :scale, scaling(o2.params.scaling, u2[:V]))
    f = KooijmanSLAPhotosynthesis();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

@testset "N assimilation" begin
    global o1, p1, u1, du1, o2, p2, u2, du2, f = nfactory();
    uptake_nitrogen(f, o1, u1)

    @testset "N assimilation depends on soil N concentration sublinealy" begin
        global a = uptake_nitrogen(f, o1, u1)
        # increase soil N
        o1.vars.assimilation.X_NO *= 2
        global b = uptake_nitrogen(f, o1, u1)
        # More N, more uptake
        @test a < b
        # But double N has less than double uptake
        @test 2a > b
    end

    @testset "N assimilation depends on root surface scaling" begin
        global a = uptake_nitrogen(f, o1, u1)
        o1.vars.scale /= 2
        global b = uptake_nitrogen(f, o1, u1)
        @test a > b
    end

    # @testset "N assimilation depends inversely on shoot surface scaling, but maybe shouldn't" begin
        # Shoot scaling inversley affects water uptake, which affects N concentration.
        # The larger the shoot is, the smaller its surface scaling coefficient. 
        # This means that except at the very start of growth, less N is assimilated by 
        # the roots the larger the shoot gets. Absolutely, not realtively. That doesn't make sense.

        # a = uptake_nitrogen(f, o1, u1)
        # o2.vars.scale /= 2
        # b = uptake_nitrogen(f, o1, u1)
        # @test a < b
    # end

    @testset "N assimilation depends on structure linearly (ignoring scaling)" begin
        global a = uptake_nitrogen(f, o1, u1)
        u1[:V] *= 2
        global b = uptake_nitrogen(f, o1, u1)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "N assimilation flux is added correctly" begin
        @test uptake_nitrogen(f, o1, u1) > 0.0u"μmol/s"
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        global c1 = sum(du1)
        global c2 = sum(du2)
        # @test n1 == o1.J[:N,:ass]
        # @test n1 == uptake_nitrogen(f, o1, u1)
    end

    @testset "N assimilation flux is merged correctly" begin

        set_rej!(o1, o2) = begin
            c_trans = 1.9oneunit(o1.J1[:1,:1])
            o1.J[:C,:tra] = o2.J[:C,:tra] = c_trans
            o1.J[:C,:rej] = o2.J[:C,:rej] = -c_trans
        end

        global o1a, p1a, u1a, du1a, o2a, p2a, u2a, du2a, fa = nfactory();
        set_rej!(o1a, o2a)
        sum_flux!(du1a, o1a, 0); sum_flux!(du2a, o2a, 0)

        global o1b, p1b, u1b, du1b, o2b, p2b, u2b, du2b, fa = nfactory();
        set_rej!(o1b, o2b)
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)

        @test o1a.J[:C,:tra] == o1b.J[:C,:tra]
        @test du1a[:C] == du1b[:C]
        @test du1a[:N] == du1b[:N]
        @test du1a[:E] == du1b[:E]
        @test sum(du2a) == sum(du2b)

        global uptake_n = uptake_nitrogen(f, o1b, u1b)
        assimilation!(f, o1b, u1b)
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)

        @test du1a[:C] > du1b[:C] 
        @test du1a[:N] < du1b[:N] 
        @test du1a[:E] < du1b[:E] 

        # Actual C and N calcs
        global c1a, c2a, c1b, c2b = sum.((du1a, du2a, du1b, du2b))
        global n1a, n2a, n1b, n2b = dot.((du1a, du2a, du1b, du2b), (n_ratios,))

        global c_loss = o1b.J1[:C,:los] + o2b.J1[:C,:los]
        global n_loss = o1b.J1[:N,:los] + o2b.J1[:N,:los]

        @test c1a > c1b
        @test n1a < n1b

        @testset "N assimilation should balance" begin
            @test -c_loss != zero(c_loss)
            @test c1b + c2b ≈ -c_loss
            @test_broken upreferred(n1b + n2b) ≈ upreferred(uptake_n) - upreferred(n_loss) 
        end
    end
end



@testset "C assimilation" begin
    global o1, p1, u1, du1, o2, p2, u2, du2, f = cfactory();
    photosynthesis(f, o1, u1)

    @testset "photosynthesis depends on air C concentration sub-linealy" begin
        global a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.X_C /= 2
        global b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "photosynthesis depends inversely on air O concentration, sub-linealy" begin
        global a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.X_O /= 2
        global b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a < b
        # But half C does not half uptake
        @test 2a > b
    end

    @testset "photosynthesis depends on Light concentration sub-linealy" begin
        global a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.J_L_F /= 2
        global b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "photosynthesis depends on shoot scaling" begin
        global a = photosynthesis(f, o1, u1)
        o1.vars.scale /= 10 
        global b = photosynthesis(f, o1, u1)
        @test a > b
    end

    # @testset "C assimilation does not depend on root scaling, but maybe it should" begin
    #     # Most stomatatal conducance models incorporate root water uptake. 
    #     # Less water means closed stomata, and less CO2. 
    #     # That is not considered in this model.

    #     a = photosynthesis(f, o1, u1)
    #     o2.vars.scale /= 2
    #     b = photosynthesis(f, o1, u1)
    #     @test a == b
    # end

    @testset "photosynthesis depends on structure linearly (ignoring scaling)" begin
        global a = photosynthesis(f, o1, u1)
        u1[:V] *= 2
        global b = photosynthesis(f, o1, u1)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "C assimilation flux is added correctly" begin
        @test photosynthesis(f, o1, u1) > 0.0u"μmol/s"
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o1, 0)
        global c1 = sum(du1)
        global c2 = sum(du2)
        @test c1 == o1.J[:C,:ass]
        @test c1 == photosynthesis(f, o1, u1)
    end

    @testset "C assimilation flux is merged correctly" begin

        global o1a, p1a, u1a, du1a, o2a, p2a, u2a, du2a, fa = cfactory();
        global n_trans = 1.9oneunit(o1a.J1[:1,:1])
        o1a.J[:N,:tra] = o2a.J[:N,:tra] = o1b.J[:N,:tra] = o2b.J[:N,:tra] = n_trans
        o1a.J[:N,:rej] = o2a.J[:N,:rej] = o1b.J[:N,:rej] = o2b.J[:N,:rej] = -n_trans
        sum_flux!(du1a, o1a, 0); sum_flux!(du2a, o2a, 0)

        global o1b, p1b, u1b, du1b, o2b, p2b, u2b, du2b, fa = cfactory();
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)

        @test du2a[:C] == du2b[:C]
        @test du2a[:N] == du2b[:N]
        @test du2a[:E] == du2b[:E]

        # assimilation only on o1b
        global uptake_c = photosynthesis(f, o1b, u1b)
        assimilation!(KooijmanSLAPhotosynthesis(), o1b, u1b)
        upreferred(o1b.J1[:N,:los])
        o1b.J
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)
        du1a
        du1b

        @test du1a[:C] < du1b[:C]
        @test du1a[:N] > du1b[:N]
        @test du1a[:E] < du1b[:E]

        global c1a, c2a, c1b, c2b = sum.((du1a, du2a, du1b, du2b))
        global n1a, n2a, n1b, n2b = dot.((du1a, du2a, du1b, du2b), (n_ratios,))
        global c_loss = o1b.J1[:C,:los] + o2b.J1[:C,:los]
        global n_loss = o1b.J1[:N,:los] + o2b.J1[:N,:los]

        @testset "C assimilation should have gained total cmols" begin
            @test c1a < c1b
            @test c2a == c2b
        end

        @testset "C assimilation should balance" begin
            @test upreferred(c1b + c2b) ≈ upreferred(-c_loss + uptake_c)
            @test upreferred(n1b + n2b) ≈ upreferred(-n_loss)
        end

    end
end
nothing
