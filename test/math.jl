using DynamicEnergyBudgets, Unitful, Test

using DynamicEnergyBudgets: rate_formula, half_saturation

@testset "Half saturation curve" begin
    mx = 10.0
    half = 5.0
    x = 5.0
    @test half_saturation(mx, half, x) == 5.0

    # Maximum
    x = 100000000000000000000000000.0
    @test half_saturation(mx, half, x) == 10.0

    # Minimum
    x = 0.0
    @test half_saturation(mx, half, x) == 0.0
end

@testset "Rate CNE" begin
    j_E_mai    = 0.003
    y_E_CH_NO  = 1.0/1.5
    y_E_EN     = 1.0/0.7
    y_V_E      = 0.7
    kap_soma   = 0.6
    r          = 0.001
    reserve    = (1.0, 1.0, 1.0)
    turnover   = ((0.2, 0.2, 0.2) .* 0.05)
    su = ParallelComplementarySU()

    rf = rate_formula(r, su, reserve, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)

    @test rf == 0.002874195250659631

    @testset "test rate with units" begin
        j_E_mai = 0.003u"mol/mol*d^-1"
        turnover = ((0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1") .* 0.05)
        mass = (1.0u"mol/mol", 1.0u"mol/mol", 1.0u"mol/mol")
        r = 0.001u"mol/mol*d^-1"
        rf = rate_formula(r, su, reserve, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)
        @test rf == 0.002874195250659631u"mol/mol*d^-1"
    end
end

@testset "Rate CN" begin
    j_E_mai    = 0.003
    y_E_CH_NO  = 1.0/1.5
    y_E_EN     = 1.0/0.7
    y_V_E      = 0.7
    kap_soma   = 0.6
    r          = 0.001
    reserve    = (10.0, 1.0)
    turnover   = ((0.2, 0.2) .* 0.05)
    su = ParallelComplementarySU()

    rf = rate_formula(r, su, reserve, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)

    # @test rf == 0.002874195250659631

    @testset "test rate with units" begin
        j_E_mai = 0.003u"mol/mol*d^-1"
        turnover = ((0.2u"mol/mol*d^-1", 0.2u"mol/mol*d^-1") .* 0.05)
        mass = (1.0u"mol/mol", 1.0u"mol/mol")
        r = 0.001u"mol/mol*d^-1"
        rf = rate_formula(r, su, reserve, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, kap_soma)
        #@test rf == 0.002874195250659631u"mol/mol*d^-1"
    end
end

nothing
