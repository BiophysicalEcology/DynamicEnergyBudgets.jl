using Revise
using Unitful
using DynamicEnergyBudgets
using DynamicEnergyBudgets: reuse_rejected!, dissipation!, translocate!, product!, 
                            maintenence!, growth!, sum_flux!, reserve_drain!, reserve_loss!,
                            maturity!, metabolism!, catabolism!, assimilation!, translocation!,
                            scaling, P, V, M, C, N, E, EE, CN, STATELEN, ass, gro, mai, rep, rej, tra, cat, rej, los,

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

sumstate(du, u) = sum(du[[P,V,M,E]]), du[C], du[N]
# sumstate(du, u) = sum(du[[1,2,5]]), du[3], du[4]

sum_c_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[C,los] + o2.J1[C,los]) * o1.shared.y_E_CH_NO
sum_n_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[N,los] + o2.J1[N,los]) * o1.shared.y_E_EN
sum_c_loss(o) = o.J1[E,los] + o.J1[C,los] * o.shared.y_E_CH_NO
sum_n_loss(o) = o.J1[E,los] + o.J1[N,los] * o.shared.y_E_EN

function factory()
    o = Organ()
    o.vars.rate = 0.1u"d^-1"
    p = o.params
    u = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    du = fill(0.0u"mol*hr^-1", 6)
    o, p, u, du
end

@testset "reserve drain" begin
    o = Organ();
    reserve_drain!(o, gro, 1.0u"mol*hr^-1", 0.4)
    @test o.J[C,gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_CH_NO
    @test o.J[N,gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_EN
    @test o.J[E,gro] ≈ 1.0u"mol*hr^-1" * 0.4 
end

@testset "growth is balanced" begin
    o, p, u, du = factory();
    o.vars.θE = 0.621

    catabolism!(o, u)
    growth!(o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "product is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, u)
    product!(o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "maintenence is balanced" begin
    o, p, u, du = factory()
    o.vars.θE = 0.621

    catabolism!(o, u)
    maintenence!(o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "maturity is balanced" begin
    o, p, u, du = factory();
    o.vars.θE = 0.621

    catabolism!(o, u)
    f = Maturity()
    maturity!(f, o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss

    # o, p, u, du = factory()
    # catabolism!(o, u)
    # u = StatePVCNE()
    # maturity!(f, o, u)
    # sum_flux!(du, 0, o)
    # m, c, n = sumstate(du, u)

    # c_loss = o.J1[E,los] + o.J1[C,los] * o.shared.y_E_CH_NO
    # n_loss = o.J1[E,los] + o.J1[N,los] * o.shared.y_E_EN
    # @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "all dissipation is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    catabolism!(o, u)
    dissipation!(o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "all metabolism is balanced" begin
    o, p, u, du = factory()
    o.vars.rate = 0.1u"d^-1"
    o.vars.θE = 0.621

    metabolism!(o, u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss

    # Works the second time?

    # reset loss because it's additive
    o.J1[E,los] = o.J1[C,los] = o.J1[N,los] = zero(o.J1[N,los])

    metabolism!(o,u)
    sum_flux!(du, o, 0)
    m, c, n = sumstate(du, u)

    c_loss = sum_c_loss(o)
    n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test m + c * o.shared.y_E_CH_NO ≈ -c_loss
    @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "translocation is balanced" begin
    o1, p1, u1, du1 = factory();
    o2, p2, u2, du2 = factory();
    o1.vars.θE = 0.621
    o2.vars.θE = 0.621
    o = o1; nothing
    on = o2; nothing
    o1.J1[E,cat] = 2oneunit(o1.J1[1,1])
    o2.J1[E,cat] = 2oneunit(o2.J1[1,1])

    translocate!(o1, o2)
    translocate!(o2, o1)
    @test maximum(abs.(o1.J)) != zero(o1.J[1,1])
    @test maximum(abs.(o2.J)) != zero(o2.J[1,1])
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    m1, c1, n1 = sumstate(du1, u1)
    m2, c2, n2 = sumstate(du2, u2)

    c_loss = sum_c_loss(o1, o2)
    n_loss = sum_n_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test m1 + m2 + (c1 + c2) * o1.shared.y_E_CH_NO ≈ -c_loss
    @test m1 + m2 + (n1 + n2) * o1.shared.y_E_EN ≈ -n_loss
end

@testset "rejection is balanced" begin
    _, _, u1, du1 = factory();
    _, _, u2, du2 = factory();
    o1 = Organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8))
    o2 = Organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8))
    p1 = o1.params; p2 = o2.params
    u2 .*= 2.7
    o1.J1[C,rej] = 2oneunit(o1.J1[1,1])
    o2.J1[C,rej] = 2oneunit(o2.J1[1,1])
    o1.J1[N,rej] = oneunit(o1.J1[1,1])
    o2.J1[N,rej] = oneunit(o2.J1[1,1])

    reuse_rejected!(o1, o2)
    reuse_rejected!(o2, o1)
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    m1, c1, n1 = sumstate(du1, u1)
    m2, c2, n2 = sumstate(du2, u2)

    c_loss = sum_c_loss(o1, o2)
    n_loss = sum_n_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test m1 + m2 + (c1 + c2) * o1.shared.y_E_CH_NO == convert(typeof(1.0u"mol/hr"), -c_loss) 
    @test m1 + m2 + (n1 + n2) * o1.shared.y_E_EN == convert(typeof(1.0u"mol/hr"), -n_loss)
end
