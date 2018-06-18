using Revise
using Unitful
using DynamicEnergyBudgets
using DynamicEnergyBudgets: reuse_rejected!, translocate!, assimilation!, translocation!, sumflux!,
                            scaling, uptake_nitrogen, photosynthesis, Kooijman_NH4_NO3Assimilation

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

sumstate(du, u::StatePVMCNE) = sum(du[[1,2,3,6]]), du[4], du[5]
sumstate(du, u::StatePVCNE) = sum(du[[1,2,5]]), du[3], du[4]

sum_C_loss(o1, o2) = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:C,:los] + o2.J1[:C,:los]) * o1.shared.y_E_CH_NO
sum_N_loss(o1, o2) = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:N,:los] + o2.J1[:N,:los]) * o1.shared.y_E_EN

function Nfactory()
    o1 = Organ(vars=Vars(assimilation=NitrogenVars()));;
    p1 = o1.params
    o2 = Organ()
    p2 = o2.params;
    u1 = StatePVMCNE(9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol")
    u2 = StatePVMCNE(9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol")
    u2 .*= 2.7
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    p1.y_EC_ECT = p1.y_EN_ENT = p2.y_EC_ECT = p2.y_EN_ENT = 0.8
    o1.vars.scale = scaling(o1.params.scaling, o1.state.V)
    o2.vars.scale = scaling(o2.params.scaling, o2.state.V)
    f = N_Assimilation();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

function Cfactory()
    o1 = Organ(vars=Vars(assimilation=CarbonVars()));;
    p1 = o1.params;
    o2 = Organ();
    p2 = o2.params;
    u1 = StatePVMCNE(9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol")
    u2 = StatePVMCNE(9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol")
    u2 .*= 1.9
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    p1.y_EC_ECT = p1.y_EN_ENT = p2.y_EC_ECT = p2.y_EN_ENT = 0.9
    o1.vars.scale = scaling(o1.params.scaling, o1.state.V)
    o2.vars.scale = scaling(o2.params.scaling, o2.state.V)
    f = KooijmanSLAPhotosynthesis();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

@testset "N assimilation" begin
    o1, p1, u1, du1, o2, p2, u2, du2, f = Nfactory();
    uptake_nitrogen(f, o1, o2)

    @testset "N assimilation depends on soil N concentration sublinealy" begin
        a = uptake_nitrogen(f, o1, o2)
        # increase soil N
        o1.vars.assimilation.X_NO *= 2
        b = uptake_nitrogen(f, o1, o2)
        # More N, more uptake
        @test a < b
        # But double N has less than double uptake
        @test 2a > b
    end

    @testset "N assimilation depends on root surface scaling" begin
        a = uptake_nitrogen(f, o1, o2)
        o1.vars.scale /= 2
        b = uptake_nitrogen(f, o1, o2)
        @test a > b
    end

    @testset "N assimilation depends inversely on shoot surface scaling, but maybe shouldn't" begin
        # Shoot scaling inversley affects water uptake, which affects N concentration.
        # The larger the shoot is, the smaller its surface scaling coefficient. 
        # This means that except at the very start of growth, less N assimilated by 
        # the roots the larger the shoot gets. That doesn't seem to make sense.

        a = uptake_nitrogen(f, o1, o2)
        o2.vars.scale /= 2
        b = uptake_nitrogen(f, o1, o2)
        @test a < b
    end

    @testset "N assimilation depends on structure linearly (ignoring scaling)" begin
        a = uptake_nitrogen(f, o1, o2)
        o1.state.V *= 2
        b = uptake_nitrogen(f, o1, o2)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "assimilation flux is added correctly" begin
        @test uptake_nitrogen(f, o1, o2) > 0.0u"μmol/s"
        uptake_nitrogen(f, o1, o2)
        assimilation!(f, o1, o2, u1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1, C1, N1 = sumstate(du1, u1)
        m2, C2, N2 = sumstate(du2, u2)
        @test N1 == o1.J[:N,:ass]
        @test N1 == uptake_nitrogen(f, o1, o2)
    end

    @testset "assimilation flux is merged correctly" begin
        o1, p1, u1, du1, o2, p2, u2, du2, f = Nfactory();
        o1.J1[:C,:rej] = 2.3oneunit(o1.J1[1,1])
        o2.J1[:C,:rej] = 2.1oneunit(o2.J1[1,1])
        o1.J1[:N,:rej] = 1.9oneunit(o1.J1[1,1])
        o2.J1[:N,:rej] = 22.3oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2)
        reuse_rejected!(o2, o1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1a, C1a, N1a = sumstate(du1, u1)
        m2a, C2a, N2a = sumstate(du2, u2)

        o1, p1, u1, du1, o2, p2, u2, du2, f = Nfactory()
        o1.J1[:C,:rej] = 2.3oneunit(o1.J1[1,1])
        o2.J1[:C,:rej] = 2.1oneunit(o2.J1[1,1])
        o1.J1[:N,:rej] = 1.9oneunit(o1.J1[1,1])
        o2.J1[:N,:rej] = 22.3oneunit(o2.J1[1,1])

        uptake_n = uptake_nitrogen(f, o1, o2)

        reuse_rejected!(o1, o2)
        reuse_rejected!(o2, o1)
        assimilation!(f, o1, o2, u1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1, C1, N1 = sumstate(du1, u1)
        m2, C2, N2 = sumstate(du2, u2)

        @test_broken C1a == C1 # Assimilation should not have added any C
        @test C2a == C2
        @test N1a != N1 # Assimilation should have added some N
        @test N2a == N2
        @test m1 != m1a # Assimilation should have added some reserve
        @test m2 == m2a

        c_loss = sum_C_loss(o1, o2)
        n_loss = sum_N_loss(o1, o2)
        @test -c_loss != zero(c_loss)
        @test upreferred(m1 + m2 + (C1 + C2) * o1.shared.y_E_CH_NO) ≈ upreferred(-c_loss)
        @test upreferred(m1 + m2 + (N1 + N2) * o1.shared.y_E_EN) ≈ -n_loss + uptake_n /o1.shared.n_N_E

    end
end



@testset "C assimilation" begin
    o1, p1, u1, du1, o2, p2, u2, du2, f = Cfactory();
    photosynthesis(f, o1, o2)

    @testset "C assimilation depends on air C concentration sub-linealy" begin
        a = photosynthesis(f, o1, o2);
        # increase soil C
        o1.vars.assimilation.X_C /= 2
        b = photosynthesis(f, o1, o2)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "C assimilation depends inversely on air O concentration, sub-linealy" begin
        a = photosynthesis(f, o1, o2);
        # increase soil C
        o1.vars.assimilation.X_O /= 2
        b = photosynthesis(f, o1, o2)
        # Less C, less uptake
        @test a < b
        # But half C does not half uptake
        @test 2a > b
    end

    @testset "C assimilation depends on Light concentration sub-linealy" begin
        a = photosynthesis(f, o1, o2);
        # increase soil C
        o1.vars.assimilation.J_L_F /= 2
        b = photosynthesis(f, o1, o2)
        # Less C, less uptake
        @test a > b
        # But half C does not half uptake
        @test a < 2b
    end

    @testset "C assimilation depends on shoot scaling" begin
        a = photosynthesis(f, o1, o2)
        o1.vars.scale /= 10 
        b = photosynthesis(f, o1, o2)
        @test a > b
    end

    @testset "C assimilation does not depend on root scaling, but maybe it should" begin
        # Most stomatatl conducance models incorporate root water uptake. 
        # Less water means closed stomata, and less CO2. 
        # That is not considered in this model.

        a = photosynthesis(f, o1, o2)
        o2.vars.scale /= 2
        b = photosynthesis(f, o1, o2)
        @test a == b
    end

    @testset "C assimilation depends on structure linearly (ignoring scaling)" begin
        a = photosynthesis(f, o1, o2)
        o1.state.V *= 2
        b = photosynthesis(f, o1, o2)
        # Lower scaling, less uptake
        @test 2a == b
    end

    @testset "C assimilation flux is added correctly" begin
        @test photosynthesis(f, o1, o2) > 0.0u"μmol/s"
        photosynthesis(f, o1, o2)
        assimilation!(f, o1, o2, u1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1, C1, N1 = sumstate(du1, u1)
        m2, C2, N2 = sumstate(du2, u2)
        @test C1 == o1.J[:C,:ass]
        @test C1 == photosynthesis(f, o1, o2)
    end

    @testset "C assimilation flux is merged correctly" begin

        o1, p1, u1, du1, o2, p2, u2, du2, f = Cfactory();
        o1.J1[:C,:rej] = 1.3oneunit(o1.J1[1,1])
        o2.J1[:C,:rej] = 3.0oneunit(o2.J1[1,1])
        o1.J1[:N,:rej] = 2.4oneunit(o1.J1[1,1])
        o2.J1[:N,:rej] = 2.9oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2)
        reuse_rejected!(o2, o1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1a, C1a, N1a = sumstate(du1, u1)
        m2a, C2a, N2a = sumstate(du2, u2)

        o1, p1, u1, du1, o2, p2, u2, du2, f = Cfactory();
        o1.J1[:C,:rej] = 1.3oneunit(o1.J1[1,1])
        o2.J1[:C,:rej] = 3.0oneunit(o2.J1[1,1])
        o1.J1[:N,:rej] = 2.4oneunit(o1.J1[1,1])
        o2.J1[:N,:rej] = 2.9oneunit(o2.J1[1,1])

        reuse_rejected!(o1, o2)
        reuse_rejected!(o2, o1)
        uptake_c = photosynthesis(f, o1, o2)
        assimilation!(f, o1, o2, u1)
        sumflux!(du1, o1, 0)
        sumflux!(du2, o2, 0)
        m1, C1, N1 = sumstate(du1, u1);
        m2, C2, N2 = sumstate(du2, u2);

        @test C1a != C1 # Assimilation should have added some C
        @test C2a == C2
        @test_broken N1a == N1 # Assimilation shouldnt have added any N
        @test N2a == N2
        @test m1 != m1a # Assimilation should have added some reserve
        @test m2 == m2a

        c_loss = sum_C_loss(o1, o2)
        n_loss = sum_N_loss(o1, o2)
        @test -c_loss != zero(c_loss)
        @test upreferred(m1 + m2 + (C1 + C2) * o1.shared.y_E_CH_NO) ≈ upreferred(-c_loss + uptake_c * o1.shared.y_E_CH_NO)
        @test m1 + m2 + (N1 + N2) * o1.shared.y_E_EN ≈ -n_loss

    end
end
