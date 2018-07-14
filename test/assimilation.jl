using Revise
using Unitful
using DynamicEnergyBudgets
using DynamicEnergyBudgets: reuse_rejected!, translocate!, assimilation!, translocation!, sum_flux!,
                            scaling, uptake_nitrogen, photosynthesis, Kooijman_NH4_NO3Assimilation,
                            P, V, M, C, N, E, ass, rej, los

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

sumstate(du, u) = sum(du[[P,V,M,E]]), du[C], du[N]

sum_c_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[C,los] + o2.J1[C,los]) * o1.shared.y_E_CH_NO
sum_n_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[N,los] + o2.J1[N,los]) * o1.shared.y_E_EN

@testset "N assimilation" begin
    o1, p1, u1, du1, o2, p2, u2, du2, f = nfactory();
    uptake_nitrogen(f, o1, u1)

    @testset "N assimilation depends on soil N concentration sublinealy" begin
        a = uptake_nitrogen(f, o1, u1)
        # increase soil N
        o1.vars.assimilation.X_NO *= 2
        b = uptake_nitrogen(f, o1, u1)
        # More N, more uptake
        @test a < b
        # But double N has less than double uptake
        @test 2a > b
    end

    @testset "N assimilation depends on root surface scaling" begin
        a = uptake_nitrogen(f, o1, u1)
        o1.vars.scale /= 2
        b = uptake_nitrogen(f, o1, u1)
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
        a = uptake_nitrogen(f, o1, u1)
        u1[V] *= 2
        b = uptake_nitrogen(f, o1, u1)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "N assimilation flux is added correctly" begin
        @test uptake_nitrogen(f, o1, u1) > 0.0u"μmol/s"
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        m1, c1, n1 = sumstate(du1, u1)
        m2, c2, n2 = sumstate(du2, u2)
        @test n1 == o1.J[N,ass]
        @test n1 == uptake_nitrogen(f, o1, u1)
    end

    @testset "N assimilation flux is merged correctly" begin

        # run without assimilation
        o1, p1, u1, du1, o2, p2, u2, du2, f = nfactory();
        o1.J1[C,rej] = 2.3oneunit(o1.J1[1,1])
        o2.J1[C,rej] = 2.1oneunit(o2.J1[1,1])
        o1.J1[N,rej] = 1.9oneunit(o1.J1[1,1])
        o2.J1[N,rej] = 22.3oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2, 1.0)
        reuse_rejected!(o2, o1, 1.0)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        m1a, c1a, n1a = sumstate(du1, u1)
        m2a, c2a, n2a = sumstate(du2, u2)

        # run with assimilation
        o1, p1, u1, du1, o2, p2, u2, du2, f = nfactory()
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
        m1, c1, n1 = sumstate(du1, u1)
        m2, c2, n2 = sumstate(du2, u2)
        
        # compare
        @testset "N assimilation should have moved some C to general, but only in o1" begin
            @test c1a > c1 
            @test c2a == c2
        end
        @testset "N assimilation should have added some N, but only in o1" begin
            @test n1a < n1 
            @test n2a == n2
        end
        @testset "N assimilation should have added some reserve, but only in o1" begin
            @test m1a < m1 
            @test m2a == m2
        end

        @testset "N assimilation should balance when converted to cmols" begin
            c_loss = sum_c_loss(o1, o2)
            n_loss = sum_n_loss(o1, o2)
            @test -c_loss != zero(c_loss)
            @test upreferred(m1 + m2 + (c1 + c2 + (uptake_n / o1.shared.n_N_N)) * o1.shared.y_E_CH_NO) ≈ upreferred(-c_loss)
            @test upreferred(m1 + m2 + (n1 + n2) * o1.shared.y_E_EN) + n_loss ≈ upreferred(uptake_n * o1.shared.y_E_EN)
        end
    end
end



@testset "C assimilation" begin
    o1, p1, u1, du1, o2, p2, u2, du2, f = cfactory();
    photosynthesis(f, o1, u1)

    @testset "C assimilation depends on air C concentration sub-linealy" begin
        a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.X_C /= 2
        b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "C assimilation depends inversely on air O concentration, sub-linealy" begin
        a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.X_O /= 2
        b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a < b
        # But half C does not half uptake
        @test 2a > b
    end

    @testset "C assimilation depends on Light concentration sub-linealy" begin
        a = photosynthesis(f, o1, u1);
        # increase soil C
        o1.vars.assimilation.J_L_F /= 2
        b = photosynthesis(f, o1, u1)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "C assimilation depends on shoot scaling" begin
        a = photosynthesis(f, o1, u1)
        o1.vars.scale /= 10 
        b = photosynthesis(f, o1, u1)
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

    @testset "C assimilation depends on structure linearly (ignoring scaling)" begin
        a = photosynthesis(f, o1, u1)
        u1[V] *= 2
        b = photosynthesis(f, o1, u1)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "C assimilation flux is added correctly" begin
        @test photosynthesis(f, o1, u1) > 0.0u"μmol/s"
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o1, 0)
        m1, c1, n1 = sumstate(du1, u1)
        m2, c2, n2 = sumstate(du2, u2)
        @test c1 == o1.J[C,ass]
        @test c1 == photosynthesis(f, o1, u1)
    end

    @testset "C assimilation flux is merged correctly" begin

        # run without assimilation
        o1, p1, u1, du1, o2, p2, u2, du2, f = cfactory();
        o1.J1[C,rej] = 1.3oneunit(o1.J1[1,1])
        o2.J1[C,rej] = 3.0oneunit(o2.J1[1,1])
        o1.J1[N,rej] = 2.4oneunit(o1.J1[1,1])
        o2.J1[N,rej] = 2.9oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2, 1.0)
        reuse_rejected!(o2, o1, 1.0)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        m1a, c1a, n1a = sumstate(du1, u1)
        m2a, c2a, n2a = sumstate(du2, u2)

        # run with assimilation
        o1, p1, u1, du1, o2, p2, u2, du2, f = cfactory();
        o1.J1[C,rej] = 1.3oneunit(o1.J1[1,1])
        o2.J1[C,rej] = 3.0oneunit(o2.J1[1,1])
        o1.J1[N,rej] = 2.4oneunit(o1.J1[1,1])
        o2.J1[N,rej] = 2.9oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2, 1.0)
        reuse_rejected!(o2, o1, 1.0)
        uptake_c = photosynthesis(f, o1, u1)
        assimilation!(f, o1, u1)
        sum_flux!(du1, o1, 0)
        sum_flux!(du2, o2, 0)
        m1, c1, n1 = sumstate(du1, u1);
        m2, c2, n2 = sumstate(du2, u2);

        # compare
        @testset "C assimilation should have added some C, but only in o1" begin
            @test c1a < c1
            @test c2a == c2
        end
        @testset "C assimilation should have moved some N to general, but only in o1" begin
            @test n1a > n1
            @test n2a == n2
        end
        @testset "C assimilation should have added some general reserve, but only in o1" begin
            @test m1a < m1 # Assimilation should have added some reserve
            @test m2a == m2
        end

        @testset "C assimilation should balance when converted to cmols" begin
            # TODO actually use cmols insteda of reserve mols
            c_loss = sum_c_loss(o1, o2)
            n_loss = sum_n_loss(o1, o2)
            @test -c_loss != zero(c_loss)
            @test upreferred(m1 + m2 + (c1 + c2) * o1.shared.y_E_CH_NO) ≈ upreferred(-c_loss + uptake_c * o1.shared.y_E_CH_NO)
            @test m1 + m2 + (n1 + n2) * o1.shared.y_E_EN ≈ -n_loss
        end

    end
end
