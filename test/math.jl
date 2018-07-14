using DynamicEnergyBudgets
using Unitful
using DynamicEnergyBudgets: catabolic_flux

@testset "Math" begin
    @test half_saturation(10.0, 5.0, 5.0) ≈ 5.0

    @test synthesizing_unit(2.0, 2.0) ≈ 1.333333333333333333
    @test synthesizing_unit(Inf, 2.0) == 2.0
    @test synthesizing_unit(0.0, 2.0) == 0.0

    (a, b, c) = stoich_merge(2.0, 2.0, 1.0, 1.0)
    @test a ≈ 0.6666666666666666666
    @test b ≈ 0.6666666666666666666
    @test c ≈ 1.3333333333333333333

    (a, b, c) = catabolic_flux((0.1,0.2,0.3), (0.2,0.3,0.4), 0.1)
    @test a ≈ 0.01
    @test b ≈ 0.04
    @test c ≈ 0.09
end

@testset "Rate" begin
    j_E_mai    = 0.003
    y_E_CH_NO  = 1.0/1.5
    y_E_EN     = 1.0/0.7
    y_V_E      = 0.7
    kap_soma   = 0.6
    r          = 0.001
    M          = (1.0, 1.0, 1.0)
    A_turnover = ((0.2, 0.2, 0.2) .* 0.05)
    rf         = rate_formula(r, M, A_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)

    @test rf == 0.002874195250659631

    @testset "test rate with units" begin
        j_E_mai   = 0.003u"mol/mol*d^-1"
        A_turnover = ((0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1") .* 0.05)
        M = (1.0u"mol/mol", 1.0u"mol/mol", 1.0u"mol/mol")
        r = 0.001u"mol/mol*d^-1"
        @test_nowarn catabolic_flux(M, A_turnover, r)
        rf = rate_formula(r, M, A_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)
        @test rf == 0.002874195250659631u"mol/mol*d^-1"
    end
end
