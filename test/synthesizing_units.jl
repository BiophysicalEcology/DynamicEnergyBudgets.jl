using DynamicEnergyBudgets, Test

@testset "ParallelComplementarySU" begin
    @test synthesizing_unit(ParallelComplementarySU(), 2.0, 2.0) ≈ 1.333333333333333333
    @test synthesizing_unit(ParallelComplementarySU(), 1000000000000000000000000.0, 2.0) == 2.0
    @test synthesizing_unit(ParallelComplementarySU(), 0.0, 2.0) == 0.0
end

@testset "MinimumRuleSU" begin
    @test synthesizing_unit(MinimumRuleSU(), 2.0, 2.0) == 2.0
    @test synthesizing_unit(MinimumRuleSU(), 1000000000000000000000000.0, 2.0) == 2.0
    @test synthesizing_unit(MinimumRuleSU(), 0.0, 2.0) == 0.0
end

@testset "KfamilySU" begin
    @test synthesizing_unit(KfamilySU(1.0), 2.0, 2.0) == 1.0
    @test synthesizing_unit(KfamilySU(2.0), 2.0, 2.0) > 1.0
    @test synthesizing_unit(KfamilySU(0.5), 2.0, 2.0) < 1.0

    @test synthesizing_unit(KfamilySU(2.0), 1000000000000000000000000.0, 2.0) ≈ 2.0
    @test synthesizing_unit(KfamilySU(1.0), 1000000000000000000000000.0, 2.0) ≈ 2.0
    @test synthesizing_unit(KfamilySU(0.5), 1000000000000000000000000.0, 2.0) ≈ 2.0

    @test synthesizing_unit(KfamilySU(2.0), 0.0, 2.0) == 0.0
    @test synthesizing_unit(KfamilySU(1.0), 0.0, 2.0) == 0.0
    @test synthesizing_unit(KfamilySU(0.5), 0.0, 2.0) == 0.0
end

@testset "SU comparisions" begin

    # Kfamiy ≈ ParallelComplementary when k ≈ 1.71. See Ledder et al, k family section. TODO: use exact page when published.
    @test synthesizing_unit(KfamilySU(1.71), 2.0, 2.0) ≈ synthesizing_unit(ParallelComplementarySU(), 2.0, 2.0) atol=0.001
    @test synthesizing_unit(KfamilySU(1.71), 3.0, 2.0) ≈ synthesizing_unit(ParallelComplementarySU(), 3.0, 2.0) atol=0.01

    # Kfamily cant handle large K falues because of floating point power limitations.
    # So we can't test the k → ∞ case where Kfamily == MinimumRule
    # @test synthesizing_unit(KfamilySU(10000.0), 2.0, 2.0) == synthesizing_unit(MinimumRuleSU(), 2.0, 2.0) == 2.0
end


@testset "Stoichiometric reserve merging" begin
    (a, b, c) = stoich_merge(ParallelComplementarySU(), 1.0, 1.0, 2.0, 2.0)

    @test a ≈ 0.3333333333333333333
    @test b ≈ 0.3333333333333333333
    @test c ≈ 1.3333333333333333333

    @test sum((a, b, c)) ≈ 1.0 + 1.0
end
