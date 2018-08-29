
sum_n_loss(o) = o.J1[E,los] + o.J1[N,los] * o.shared.y_E_EN
sum_n_loss(o1, o2) = o1.J1[E,los] + o2.J1[E,los] + (o1.J1[N,los] + o2.J1[N,los]) * o1.shared.y_E_EN

@testset "stoich merge is balanced" begin
    @test sum(stoich_merge(1.0, 2.0, 0.7, 0.4)) == 3.0
end


@testset "reserve drain works" begin
    global o, p, u, du = factory();
    reserve_drain!(o, gro, 1.0u"mol*hr^-1", 0.4)
    @test o.J[C,gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_CH_NO
    @test o.J[N,gro] ≈ 1.0u"mol*hr^-1" * (1 - 0.4)/o.shared.y_E_EN
    @test o.J[E,gro] ≈ 1.0u"mol*hr^-1" * 0.4
end

@testset "catabolism works" begin
    global o, p, u, du = factory();
    @test o.J1[C,cat]  == zero(o.J1[1,1])
    @test o.J1[N,cat]  == zero(o.J1[1,1])
    @test o.J1[C,rej]  == zero(o.J1[1,1])
    @test o.J1[N,rej]  == zero(o.J1[1,1])
    @test o.J1[E,cat]  == zero(o.J1[1,1])
    @test o.J1[EE,cat] == zero(o.J1[1,1]) 
    @test o.J1[CN,cat] == zero(o.J1[1,1])
    @test o.J1[C,los] == zero(o.J1[1,1])
    @test o.J1[N,los] == zero(o.J1[1,1])

    catabolism!(o, u)

    @test o.J1[C,cat]  > zero(o.J1[1,1])
    @test o.J1[N,cat]  > zero(o.J1[1,1])
    @test o.J1[C,rej]  > zero(o.J1[1,1])
    @test o.J1[N,rej]  > zero(o.J1[1,1])
    @test o.J1[E,cat]  > zero(o.J1[1,1])
    @test o.J1[EE,cat] > zero(o.J1[1,1])
    @test o.J1[CN,cat] > zero(o.J1[1,1])
    @test o.J1[EE,cat] + o.J1[CN,cat] == o.J1[E,cat]
    @test o.J1[C,los] > zero(o.J1[1,1])
    @test_broken o.J1[N,los] > zero(o.J1[1,1])
end

@testset "growth is balanced" begin
    global o, p, u, du = factory();
    set_var!(o.vars, :θE, 0.621)

    catabolism!(o, u)
    global cat_loss = sum(o.J1[:, los])
    growth!(o, u)
    sum_flux!(du, o, 0)

    global c = sum(du)
    global c_loss = sum(o.J1[:,los]) - cat_loss
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "product is balanced" begin
    global o, p, u, du = factory()
    set_var!(o.vars, :θE, 0.621)

    catabolism!(o, u)
    global cat_loss = sum(o.J1[:, los])

    product!(o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)

    global c_loss = sum(o.J1[:,los]) - cat_loss
    # n_loss = sum_n_loss(o)
    @test c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end



@testset "maintenence is balanced" begin
    global o, p, u, du = factory()
    set_var!(o.vars, :θE, 0.621)

    catabolism!(o, u)
    global cat_loss = sum(o.J1[:, los])

    maintenence!(o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)

    global c_loss = sum(o.J1[:,los]) - cat_loss
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "maturity is balanced" begin
    global o, p, u, du = factory();
    set_var!(o.vars, :θE, 0.621)

    catabolism!(o, u)
    global cat_loss = sum(o.J1[:, los])

    global f = Maturity()
    maturity!(f, o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)

    global c_loss = sum(o.J1[:,los]) - cat_loss
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "all dissipation is balanced" begin
    global o, p, u, du = factory()
    set_var!(o.vars, :rate, 0.1u"d^-1")
    set_var!(o.vars, :θE, 0.621)

    catabolism!(o, u)
    global cat_loss = sum(o.J1[:, los])

    dissipation!(o, u)
    sum_flux!(du, o, 0)

    global c = sum(du)
    global c_loss = sum(o.J1[:,los]) - cat_loss
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "all metabolism is balanced" begin
    global o, p, u, du = factory()
    set_var!(o.vars, :rate, 0.1u"d^-1")
    set_var!(o.vars, :θE, 0.621)

    metabolism!(o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)

    global c_loss = sum(o.J1[:,los])
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss

    # Works the second time?

    # reset loss because it's additive
    o.J1[E,los] = o.J1[C,los] = o.J1[N,los] = zero(o.J1[N,los])

    metabolism!(o,u)
    sum_flux!(du, o, 0)
    global c = sum(du)

    global c_loss = sum(o.J1[:,los])
    # n_loss = sum_n_loss(o)
    @test -c_loss != zero(c_loss)
    @test_broken c ≈ -c_loss
    # @test m + n * o.shared.y_E_EN ≈ -n_loss
end

@testset "translocation is balanced" begin
    global o1, p1, u1, du1 = factory();
    global o2, p2, u2, du2 = factory();
    set_var!(o1.vars, :θE, 0.621)
    set_var!(o2.vars, :θE, 0.621)
    global o = o1; nothing
    global on = o2; nothing
    o1.J1[E,cat] = 2oneunit(o1.J1[1,1])
    o2.J1[E,cat] = 2oneunit(o2.J1[1,1])

    translocate!(o1, o2, 1.0)
    translocate!(o2, o1, 1.0)
    @test maximum(abs.(o1.J)) != zero(o1.J[1,1])
    @test maximum(abs.(o2.J)) != zero(o2.J[1,1])
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    global c1 = sum(du1)
    global c2 = sum(du2)

    global c_loss = sum(o1.J1[:,los]) + sum(o2.J1[:,los])
    # n_loss = sum_n_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test_broken c1 + c2 ≈ -c_loss
    # @test m1 + m2 + (n1 + n2) * o1.shared.y_E_EN ≈ -n_loss
end

@testset "rejection is balanced" begin
    global _, _, u1, du1 = factory();
    global _, _, u2, du2 = factory();
    global o1 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8))
    global o2 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8))
    global p1 = o1.params; 
    global p2 = o2.params
    u2 .*= 2.7
    o1.J1[C,rej] = 2oneunit(o1.J1[1,1])
    o2.J1[C,rej] = 2oneunit(o2.J1[1,1])
    o1.J1[N,rej] = oneunit(o1.J1[1,1])
    o2.J1[N,rej] = oneunit(o2.J1[1,1])

    reuse_rejected!(o1, o2, 1.0)
    reuse_rejected!(o2, o1, 1.0)
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    global c1 = sum(du1)
    global c2 = sum(du2)

    global c_loss = sum(o1.J1[:,los]) + sum(o2.J1[:,los])
    # n_loss = sum_n_loss(o1, o2)
    @test -c_loss != zero(c_loss)
    @test_broken c1 + c2 == -c_loss
    # @test m1 + m2 + (n1 + n2) * o1.shared.y_E_EN == convert(typeof(1.0u"mol/hr"), -n_loss)
end
nothing
