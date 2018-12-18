
@testset "stoich merge is balanced if yeilds are 2.0" begin
    @test sum(stoich_merge(2.0, 2.0, 2.0, 2.0)) ≈ 4.0
    @test sum(stoich_merge(1.0, 3.0, 2.0, 2.0)) ≈ 4.0
    @test sum(stoich_merge(100.0, 0.0, 2.0, 2.0)) ≈ 100.0
end

@testset "reserve drain and loss match" begin
    global o, p, u, du = factory();
    set_var!(o.vars, :θE, 0.4)
    reserve_drain!(o, :gro, 1.0u"mol*hr^-1")
    @test o.J[:C,:gro] ≈ -1.0u"mol*hr^-1" * (1 - 0.4)/o.params.y_E_CH_NO
    @test o.J[:N,:gro] ≈ -1.0u"mol*hr^-1" * (1 - 0.4)/o.params.y_E_EN
    @test o.J[:E,:gro] ≈ -1.0u"mol*hr^-1" * 0.4
    reserve_loss!(o, 1.0u"mol*hr^-1")
    @test o.J1[:C,:los] == -sum(o.J[:,:gro])
    @test o.J1[:N,:los] == -o.J[:,:gro] ⋅ n_ratios
end

@testset "catabolism works" begin
    global o, p, u, du = factory();
    @test o.J1[:C,:ctb]  == zero(o.J1[:1,:1])
    @test o.J1[:N,:ctb] == zero(o.J1[:1,:1])
    @test o.J1[:C,:rej]  == zero(o.J1[:1,:1])
    @test o.J1[:N,:rej]  == zero(o.J1[:1,:1])
    @test o.J1[:E,:ctb]  == zero(o.J1[:1,:1])
    @test o.J1[:EE,:ctb] == zero(o.J1[:1,:1])
    @test o.J1[:CN,:ctb] == zero(o.J1[:1,:1])
    @test o.J1[:C,:los]  == zero(o.J1[:1,:1])
    @test o.J1[:N,:los]  == zero(o.J1[:1,:1])

    catabolism!(o, u)

    @test o.J1[:C,:ctb]  > zero(o.J1[:1,:1])
    @test o.J1[:N,:ctb]  > zero(o.J1[:1,:1])
    @test o.J1[:C,:rej]  > zero(o.J1[:1,:1])
    @test o.J1[:N,:rej]  > zero(o.J1[:1,:1])
    @test o.J1[:E,:ctb]  > zero(o.J1[:1,:1])
    @test o.J1[:EE,:ctb] > zero(o.J1[:1,:1])
    @test o.J1[:CN,:ctb] > zero(o.J1[:1,:1])
    @test o.J1[:EE,:ctb] + o.J1[:CN,:ctb] == o.J1[:E,:ctb]
end

@testset "growth is balanced" begin
    global o, p, u, du = factory();
    set_scaling!(o, u)
    catabolism!(o, u)
    @test rate(o.vars) > zero(rate(o.vars)) 

    growth!(o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)
    global n = du ⋅ n_ratios

    @test c < zero(c)
    @test c ≈ -o.J1[:C,:los]
    @test n ≈ -o.J1[:N,:los]

end

@testset "maintenence is balanced" begin
    global o, p, u, du = factory()
    set_scaling!(o, u)
    catabolism!(o, u)

    maintenence!(o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)
    global n = du ⋅ n_ratios
    @test c < zero(c)
    @test n < zero(n)
    @test c ≈ -o.J1[:C,:los]
    @test n ≈ -o.J1[:N,:los]

end

@testset "maturity is balanced" begin
    global o, p, u, du = factory();
    set_scaling!(o, u)
    catabolism!(o, u)

    maturity!(Maturity(), o, u)
    sum_flux!(du, o, 0)
    global c = sum(du)
    global n = du ⋅ n_ratios
    @test c < zero(c)
    @test -c ≈ o.J1[:C,:los]
    @test -n ≈ o.J1[:N,:los]

end

@testset "all metabolism is balanced" begin
    global o, p, u, du = factory()
    metabolism!(o, u)
    sum_flux!(du, o, 0)

    global c = sum(du)
    global n = du ⋅ n_ratios
    @test c < zero(c)
    @test -c ≈ o.J1[:C,:los]
    @test -n ≈ o.J1[:N,:los]

end

@testset "translocation is balanced" begin
    global o1, p1, u1, du1 = factory();
    global o2, p2, u2, du2 = factory();
    set_var!(o1.vars, :θE, 0.621)
    set_var!(o2.vars, :θE, 0.621)
    global o = o1; nothing
    global on = o2; nothing
    o1.J1[:E,:ctb] = 2oneunit(o1.J1[:1,:1])
    o2.J1[:E,:ctb] = 2oneunit(o2.J1[:1,:1])

    translocate!(o1, o2, 1.0)
    translocate!(o2, o1, 1.0)
    @test maximum(abs.(o1.J)) != zero(o1.J[:1,:1])
    @test maximum(abs.(o2.J)) != zero(o2.J[:1,:1])
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    global c1 = sum(du1)
    global c2 = sum(du2)
    global n1 = du1 ⋅ n_ratios
    global n2 = du1 ⋅ n_ratios
    @test c1 != zero(c1)
    @test -c1 - c2 ≈ o1.J1[:C,:los] + o2.J1[:C,:los]
    @test -n1 - n2 ≈ o1.J1[:N,:los] + o2.J1[:N,:los]

end

@testset "rejection is balanced" begin
    global o1, p1, u1, du1 = factory();
    global o2, p2, u2, du2 = factory();

    # These are set to 1.0 by default which makes this test meaningless
    global o1 = construct_organ(params=Params(y_EC_ECT = 0.7, y_EN_ENT = 0.9))
    global o2 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.5))
    global p1 = o1.params; p2 = o2.params
    global sh = o1.shared
    u2 .*= 2.7
    o1.J1[:C,:rej] = 2oneunit(o1.J1[:1,:1])
    o2.J1[:C,:rej] = 2oneunit(o2.J1[:1,:1])
    o1.J1[:N,:rej] = oneunit(o1.J1[:1,:1])
    o2.J1[:N,:rej] = oneunit(o2.J1[:1,:1])

    reuse_rejected!(o1, o2, 1.0)
    reuse_rejected!(o2, o1, 1.0)
    sum_flux!(du1, o1, 0)
    sum_flux!(du2, o2, 0)
    global c1 = sum(du1)
    global c2 = sum(du2)
    global n1 = du1 ⋅ n_ratios
    global n2 = du2 ⋅ n_ratios

    @test c1 != zero(c1)
    @test -c1 - c2 ≈ o1.J1[:C,:los] + o2.J1[:C,:los]
    @test -n1 - n2 ≈ o1.J1[:N,:los] + o2.J1[:N,:los]

end

nothing
