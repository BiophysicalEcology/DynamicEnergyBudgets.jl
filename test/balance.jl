using Revise
using DynamicEnergyBudgets
using DynamicEnergyBudgets: reuse_rejected!, dissipation!, translocate!, product!, 
                            maintenence!, growth!, sumflux!, reserve_drain!, reserve_loss!,
                            maturity!, metabolism!, catabolism!, assimilation!, translocation!,
                            scaling
using Unitful

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

sumstate(du, u::StatePVMCNE) = sum(du[[1,2,3,6]]), du[4], du[5]
sumstate(du, u::StatePVCNE) = sum(du[[1,2,5]]), du[3], du[4]

sum_C_loss(o1, o2) = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:C,:los] + o2.J1[:C,:los]) * o1.shared.y_E_CH_NO
sum_N_loss(o1, o2) = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:N,:los] + o2.J1[:N,:los]) * o1.shared.y_E_EN
sum_C_loss(o) = o.J1[:E,:los] + o.J1[:C,:los] * o.shared.y_E_CH_NO
sum_N_loss(o) = o.J1[:E,:los] + o.J1[:N,:los] * o.shared.y_E_EN

function factory()
    o = Organ()
    p = o.params
    u = StatePVMCNE(9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol")
    du = fill(0.0u"mol*hr^-1", 6)
    o, p, u, du
end

@testset "reserve drain" begin
    o = Organ();
    reserve_drain!(o, :gro, 1.0u"mol*hr^-1", 0.4)
    @test o.J[:C,:gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_CH_NO
    @test o.J[:N,:gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_EN
    @test o.J[:E,:gro] ≈ 1.0u"mol*hr^-1" * 0.4 
end

@testset "growth is balanced" begin
    o, p, u, du = factory();
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    growth!(o, u)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "product is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    product!(o, u)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "maintenence is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    maintenence!(o, u)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "maturity is balanced" begin
    o, p, u, du = factory();
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    f = Maturity()
    maturity!(f, o, u)
    sumflux!(du, o, 0)
    du
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss

    # o, p, u, du = factory()
    # catabolism!(o, o.state, 1)
    # u = StatePVCNE()
    # maturity!(f, o, u)
    # sumflux!(du, 0, o)
    # m, C, N = sumstate(du, u)

    # c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.shared.y_E_CH_NO
    # n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.shared.y_E_EN
    # @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    # @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "all dissipation is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    dissipation!(o, u)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "all metabolism is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    metabolism!(o, 1)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss

    # Works the second time?

    # reset loss because it's additive
    o.J1[:E,:los] = o.J1[:C,:los] = o.J1[:N,:los] = zero(o.J1[:N,:los])

    metabolism!(o, 1)
    sumflux!(du, o, 0)
    m, C, N = sumstate(du, u)

    c_loss = sum_C_loss(o)
    n_loss = sum_N_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + C * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + N * o.shared.y_E_EN ≈ -n_loss
end

@testset "translocation is balanced" begin
    o1, p1, u1, du1 = factory();
    o2, p2, u2, du2 = factory();
    o1.vars.θE = 0.621
    o2.vars.θE = 0.621
    o = o1; nothing
    on = o2; nothing
    o1.J1[:E,:cat] = 2oneunit(o1.J1[1,1])
    o2.J1[:E,:cat] = 2oneunit(o2.J1[1,1])

    translocate!(o1, o2)
    translocate!(o2, o1)
    @test maximum(abs.(o1.J)) != zero(o1.J[1,1])
    @test maximum(abs.(o2.J)) != zero(o2.J[1,1])
    sumflux!(du1, o1, 0)
    sumflux!(du2, o2, 0)
    m1, C1, N1 = sumstate(du1, u1)
    m2, C2, N2 = sumstate(du2, u2)

    c_loss = sum_C_loss(o1, o2)
    n_loss = sum_N_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test m1 + m2 + (C1 + C2) * o1.shared.y_E_CH_NO ≈ -c_loss
    @test m1 + m2 + (N1 + N2) * o1.shared.y_E_EN ≈ -n_loss
end

@testset "rejection is balanced" begin
    o1, p1, u1, du1 = factory();
    o2, p2, u2, du2 = factory();
    u2 .*= 2.7
    o1.J1[:C,:rej] = 2oneunit(o1.J1[1,1])
    o2.J1[:C,:rej] = 2oneunit(o2.J1[1,1])
    o1.J1[:N,:rej] = oneunit(o1.J1[1,1])
    o2.J1[:N,:rej] = oneunit(o2.J1[1,1])

    p1.y_EC_ECT = p1.y_EN_ENT = p2.y_EC_ECT = p2.y_EN_ENT = 0.8
    reuse_rejected!(o1, o2)
    reuse_rejected!(o2, o1)
    sumflux!(du1, o1, 0)
    sumflux!(du2, o2, 0)
    m1, C1, N1 = sumstate(du1, u1)
    m2, C2, N2 = sumstate(du2, u2)

    c_loss = sum_C_loss(o1, o2)
    n_loss = sum_N_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test m1 + m2 + (C1 + C2) * o1.shared.y_E_CH_NO == -c_loss
    @test m1 + m2 + (N1 + N2) * o1.shared.y_E_EN == -n_loss
end
