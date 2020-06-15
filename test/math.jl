using DynamicEnergyBudgets, Unitful, Test
using DynamicEnergyBudgets: non_growth_flux

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

@testset "Flux" begin
    (a, b, c) = non_growth_flux.((0.1,0.2,0.3), (0.2,0.3,0.4), 0.1)
    @test a ≈ 0.01
    @test b ≈ 0.04
    @test c ≈ 0.09
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

    # rate_formula(r, su, rel_reserve::NTuple{3}, turnover::NTuple{3}, j_E_mai, y_E_Ea, y_E_Eb, y_V_E, κsoma) = begin
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

nothing
