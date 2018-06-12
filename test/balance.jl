using DynamicEnergyBudgets: reuse_rejected!, dissipation!, translocate!, product!, 
                            maintenence!, growth!, sum_flux!, reserve_drain!, reserve_loss!,
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


function factory()
    o = Organ()
    p = o.params
    u = StatePVMCNE()
    du = fill(0.0u"mol*hr^-1", 6)
    o, p, u, du
end

@testset "reserve drain" begin
    o = Organ()
    reserve_drain!(o, :gro, 1.0u"mol*hr^-1", 0.4)
    @test o.J[:C,:gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.params.y_E_CH_NO
    @test o.J[:N,:gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.params.y_E_EN
    @test o.J[:E,:gro] ≈ 1.0u"mol*hr^-1" * 0.4 
end

@testset "growth is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    growth!(o, u)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "product is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    product!(o, u)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "maintenence is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    maintenence!(o, u)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "maturity is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    f = Maturity()
    maturity!(f, o, u)
    sum_flux!(du, 0, o)
    du
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss

    # o, p, u, du = factory()
    # catabolism!(o, o.state, 1)
    # u = StatePVCNE()
    # maturity!(f, o, u)
    # sum_flux!(du, 0, o)
    # m, C, N = sumstate(du, u)

    # c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    # n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    # @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    # @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "all dissipation is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    catabolism!(o, o.state, 1)
    dissipation!(o, u)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "all metabolism is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    metabolism!(o, 1)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss

    # Works the second time?

    # reset loss because it's additive
    o.J1[:E,:los] = o.J1[:C,:los] = o.J1[:N,:los] = zero(o.J1[:N,:los])

    metabolism!(o, 1)
    sum_flux!(du, 0, o)
    m, C, N = sumstate(du, u)

    c_loss = o.J1[:E,:los] + o.J1[:C,:los] * o.params.y_E_CH_NO
    n_loss = o.J1[:E,:los] + o.J1[:N,:los] * o.params.y_E_EN
    @test m + C * o.params.y_E_CH_NO ≈ -c_loss
    @test m + N * o.params.y_E_EN ≈ -n_loss
end

@testset "translocation is balanced" begin
    o1, p1, u1, du1 = factory()
    o2, p2, u2, du2 = factory()
    o1.vars.θE = 0.621
    o2.vars.θE = 0.621
    o1.J1[:C,:rej], 2oneunit(o1.J1[1,1])
    o2.J1[:C,:rej], 2oneunit(o2.J1[1,1])
    o1.J1[:N,:rej], oneunit(o1.J1[1,1])
    o2.J1[:N,:rej], oneunit(o2.J1[1,1])

    translocate!(o1, o2)
    sum_flux!(du1, 0, o1)
    sum_flux!(du2, 0, o2)
    m1, C1, N1 = sumstate(du1, u1)
    m2, C2, N2 = sumstate(du2, u2)

    c_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:C,:los] + o2.J1[:C,:los]) * o1.params.y_E_CH_NO
    n_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:N,:los] + o2.J1[:N,:los]) * o1.params.y_E_EN
    @test m1 + m2 + (C1 + C2) * o1.params.y_E_CH_NO == -c_loss
    @test m1 + m2 + (N1 + N2) * o1.params.y_E_EN == -n_loss
end

@testset "rejection is balanced" begin
    o1, p1, u1, du1 = factory()
    o2, p2, u2, du2 = factory()
    o1.vars.θE = o2.vars.θE = 0.621
    o1.J1[:C,:rej] = 2oneunit(o1.J1[1,1])
    o2.J1[:C,:rej] = 2oneunit(o2.J1[1,1])
    o1.J1[:N,:rej] = oneunit(o1.J1[1,1])
    o2.J1[:N,:rej] = oneunit(o2.J1[1,1])

    p1.y_EC_ECT = p1.y_EN_ENT = p2.y_EC_ECT = p2.y_EN_ENT = 0.8
    reuse_rejected!(o1, o2)
    sum_flux!(du1, 0, o1)
    sum_flux!(du2, 0, o2)
    m1, C1, N1 = sumstate(du1, u1)
    m2, C2, N2 = sumstate(du2, u2)

    c_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:C,:los] + o2.J1[:C,:los]) * o1.params.y_E_CH_NO
    n_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:N,:los] + o2.J1[:N,:los]) * o1.params.y_E_EN
    @test m1 + m2 + (C1 + C2) * o1.params.y_E_CH_NO == -c_loss
    @test m1 + m2 + (N1 + N2) * o1.params.y_E_EN == -n_loss
end

@testset "assimilation is balanced" begin
    o, p, u, du = factory()
    o2, p2, u2, du2 = factory()
    o1.vars.θE = 0.621
    o2.vars.θE = 0.621
    o1.J1[:C,:rej], 2oneunit(o1.J1[1,1])
    o2.J1[:C,:rej], 2oneunit(o2.J1[1,1])
    o1.J1[:N,:rej], oneunit(o1.J1[1,1])
    o2.J1[:N,:rej], oneunit(o2.J1[1,1])

    u.V = oneunit(u.V)
    o.vars.scale = scaling(o.params.scaling, o.state.V)
    translocation!(o, o2)
    f = Kooijman_NH4_NO3_Assimilation()

    # DynamicEnergyBudgets.uptake_nitrogen(f::Kooijman_NH4_NO3_Assimilation, o, on) = (0.0u"mol*mol^-1*s^-1", 0.0u"mol*mol^-1*s^-1", 0.0u"mol*mol^-1*s^-1")
    
    assimilation!(o, o2, f, u)
    sum_flux!(du, 0, o1)
    sum_flux!(du2, 0, o2)
    m, C, N = sumstate(du, u)
    m2, C2, N2 = sumstate(du2, u2)

    c_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:C,:los] + o2.J1[:C,:los]) * o1.params.y_E_CH_NO
    n_loss = o1.J1[:E,:los] + o2.J1[:E,:los] + (o1.J1[:N,:los] + o2.J1[:N,:los]) * o1.params.y_E_EN
    @test m1 + m2 + (C1 + C2) * o1.params.y_E_CH_NO == -c_loss + o1.J[:E,:ass]
    @test m1 + m2 + (N1 + N2) * o1.params.y_E_EN == -(n_loss + o.J[:E,:ass] + o.J[:N,:ass] * o1.params.y_E_EN)
end

