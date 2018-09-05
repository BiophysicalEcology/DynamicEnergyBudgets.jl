using DynamicEnergyBudgets: reuse_rejected!, translocate!, assimilation!, translocation!, sum_flux!,
                            scaling, uptake_nitrogen, photosynthesis,
                            P, V, M, C, N, E, ass, rej, los

sum_n_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[N,los] + o2.J1[N,los]) * o1.shared.y_E_EN

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
        u1[V] *= 2
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
        # @test n1 == o1.J[N,ass]
        # @test n1 == uptake_nitrogen(f, o1, u1)
    end

    @testset "N assimilation flux is merged correctly" begin

        # run without assimilation
        global o1, p1, u1, du1a, o2, p2, u2, du2a, f = nfactory();
        o1.J1[C,rej] = 2.3oneunit(o1.J1[1,1])
        o2.J1[C,rej] = 2.1oneunit(o2.J1[1,1])
        o1.J1[N,rej] = 1.9oneunit(o1.J1[1,1])
        o2.J1[N,rej] = 22.3oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2, 1.0)
        reuse_rejected!(o2, o1, 1.0)
        sum_flux!(du1a, o1, 0)
        sum_flux!(du2a, o2, 0)
        global c1a = sum(du1a)
        global c2a = sum(du2a)

        # run with assimilation
        global o1, p1, u1, du1, o2, p2, u2, du2, f = nfactory()
        o1.J1[C,rej] = 2.3oneunit(o1.J1[1,1])
        o2.J1[C,rej] = 2.1oneunit(o2.J1[1,1])
        o1.J1[N,rej] = 1.9oneunit(o1.J1[1,1])
        o2.J1[N,rej] = 22.3oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2, 1.0)
        reuse_rejected!(o2, o1, 1.0)
        uptake_n = uptake_nitrogen(f, o1, u1)
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        global c1 = sum(du1)
        global c2 = sum(du2)
        
        # compare
        @testset "N assimilation should have lost C and moved some to general" begin
            @test du1a[C] > du1[C] 
            @test du2a[C] == du2[C]
        end
        @testset "N assimilation should have added N" begin
            @test du1a[N] < du1[N] 
            @test du2a[N] == du2[N] 
        end
        @testset "N assimilation should have reduced total cmols" begin
            @test_broken c1a > c1 
            @test c2a == c2
        end

        @testset "N assimilation should balance" begin
            global c_loss = sum(o1.J1[:,los]) + sum(o2.J1[:,los])
            global n_loss = sum_n_loss(o1, o2)
            @test -c_loss != zero(c_loss)
            @test_broken c1 + c2 + (uptake_n / o1.shared.n_N_N) ≈ -c_loss
            @test_broken upreferred(n1 + n2 * o1.shared.y_E_EN) + n_loss ≈ upreferred(uptake_n * o1.shared.y_E_EN)
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
        u1[V] *= 2
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
        @test c1 == o1.J[C,ass]
        @test c1 == photosynthesis(f, o1, u1)
    end

    @testset "C assimilation flux is merged correctly" begin

        global o1a, p1a, u1a, du1a, o2a, p2a, u2a, du2a, fa = cfactory();
        n_trans = 1.9oneunit(o1a.J1[1,1])
        o1a.J[N,tra] = o2a.J[N,tra] = n_trans
        o1a.J[N,rej] = o2a.J[N,rej] = -n_trans
        sum_flux!(du1a, o1a, 0); sum_flux!(du2a, o2a, 0)

        global o1b, p1b, u1b, du1b, o2b, p2b, u2b, du2b, fa = cfactory();
        n_trans = 1.9oneunit(o1a.J1[1,1])
        o1b.J[N,tra] = o2b.J[N,tra] = n_trans
        o1b.J[N,rej] = o2b.J[N,rej] = -n_trans
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)

        @test du2a[C] == du2b[C]
        @test du2a[N] == du2b[N]
        @test du2a[E] == du2b[E]

        # assimilation only on o1b
        global uptake_c = photosynthesis(f, o1b, u1b)
        assimilation!(KooijmanSLAPhotosynthesis(), o1b, u1b)
        # (o1b.J[C,ass], o1b.J[N,tra], o1b.J[E,ass], lossC, lossN) =
            # stoich_merge(uptake_c, o1b.J[N,tra], o1b.shared.y_E_CH_NO, o1b.shared.y_E_EN)
        # o1b.J1[C,los] += lossC 
        # o1b.J[N,los] += lossN
        sum_flux!(du1b, o1b, 0); sum_flux!(du2b, o2b, 0)

        # compare
        @testset "C assimilation should have added some C" begin
            @test du1a[C] < du1b[C]
        end
        @testset "C assimilation should have moved some N to general" begin
            @test du1a[N] > du1b[N]
            @test du1a[E] < du1b[E]
        end

        global cmol1a, cmol2a, cmol1, cmol2 = sum(du1a), sum(du2a), sum(du1b), sum(du2b)

        @testset "C assimilation should have gained total cmols" begin
            @test cmol1a < cmol1
            @test cmol2a == cmol2
        end

        @testset "C assimilation should balance" begin
            upreferred(uptake_c)
            upreferred(cmol1 + cmol2)
            @test_broken upreferred(cmol1 + cmol2) ≈ upreferred(uptake_c)
        end

    end
end
nothing
