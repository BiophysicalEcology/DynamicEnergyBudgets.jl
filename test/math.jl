using DynamicEnergyBudgets
using Unitful
using DynamicEnergyBudgets: catabolic_flux

@testset "Math" begin
    @test half_saturation(10.0, 5.0, 5.0) ≈ 5.0

    @test synthesizing_unit(2.0, 2.0) ≈ 1.333333333333333333
    @test synthesizing_unit(Inf, 2.0) == 2.0
    @test synthesizing_unit(0.0, 2.0) == 0.0

    global (a, b, c) = stoich_merge(1.0, 1.0, 2.0, 2.0)
    @test a ≈ 0.3333333333333333333
    @test b ≈ 0.3333333333333333333
    @test c ≈ 1.3333333333333333333
    @test sum((a, b, c)) ≈ 1.0 + 1.0

    global (a, b, c) = catabolic_flux.((0.1,0.2,0.3), (0.2,0.3,0.4), 0.1)
    @test a ≈ 0.01
    @test b ≈ 0.04
    @test c ≈ 0.09
end

@testset "Rate" begin
    global j_E_mai    = 0.003
    global y_E_CH_NO  = 1.0/1.5
    global y_E_EN     = 1.0/0.7
    global y_V_E      = 0.7
    global kap_soma   = 0.6
    global r          = 0.001
    global mass       = (1.0, 1.0, 1.0)
    global A_turnover = ((0.2, 0.2, 0.2) .* 0.05)
    global rf         = rate_formula(r, mass, A_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)

    @test rf == 0.002874195250659631

    @testset "test rate with units" begin
        global j_E_mai   = 0.003u"mol/mol*d^-1"
        global A_turnover = ((0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1") .* 0.05)
        global mass = (1.0u"mol/mol", 1.0u"mol/mol", 1.0u"mol/mol")
        global r = 0.001u"mol/mol*d^-1"
        @test_nowarn catabolic_flux.(mass, A_turnover, r)
        global rf = rate_formula(r, mass, A_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)
        @test rf == 0.002874195250659631u"mol/mol*d^-1"
    end

end

nothing
