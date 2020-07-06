using Test, DynamicEnergyBudgets, Unitful, OrdinaryDiffEq, DimensionalData
using Microclimate, Photosynthesis

u0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]u"mol" # Initial value
du0 = DimArray(fill(0.0u"mol/hr", 6), (X(Val((:VS,:CS,:NS,:VR,:CR,:NR))),))

env = MicroclimControl()
maintenance_cost = -0.001u"hr^-1"
shared=SharedParams(
    resorption_pars=nothing,
    tempcorr_pars=nothing,
    su_pars = ParallelComplementarySU(),
    core_pars = DEBCore(
        j_E_mai = maintenance_cost,
        y_V_E = 1.0,
        y_E_C = 1.0,
        y_E_N = 100.0,
        n_N_V = 0.01,
        n_N_E = 0.01,
    ),
    catabolism_pars = CatabolismCNshared(
        k = 0.5u"d^-1",
    ),
)

@testset "balance with growth and maintenance is correct" begin
    plant1 = Plant(time = 0u"hr":1u"hr":8001u"hr", 
        params=(
            Params(
                assimilation_pars = ConstantCAssim(
                    c_uptake = 0.0u"mol*mol^-1*hr^-1",
                ),
            ),  
            Params(
                assimilation_pars = ConstantNAssim(
                    n_uptake = 0.0u"mol*mol^-1*hr^-1",
                ),
            ),
        ),
        shared=shared,
        environment=env
    )

    du1 = deepcopy(du0)
    u1 = deepcopy(u0)
    plant1(du1, u1, nothing, 1u"hr")    

    @test du1[:VS] > 0.0u"mol/hr"
    @test du1[:CS] < 0.0u"mol/hr"
    @test du1[:NS] < 0.0u"mol/hr"
    @test du1[:VR] > 0.0u"mol/hr"
    @test du1[:CR] < 0.0u"mol/hr"
    @test du1[:NR] < 0.0u"mol/hr"

    # Carbon is balanced, with losses to maintenance
    @test du1[:VS] == -du1[:CS] - maintenance_cost * 1u"mol"
    @test du1[:VR] == -du1[:CR] - maintenance_cost * 1u"mol"
    @test du1[:VS] ≈ -100du1[:NS] - maintenance_cost * 1u"mol"
    @test du1[:VR] ≈ -100du1[:NR] - maintenance_cost * 1u"mol"
end

@testset "balance with assimilation is correct" begin
    plant = Plant(time = 0u"hr":1u"hr":8001u"hr", 
        params=(
            Params(
                assimilation_pars = ConstantCAssim(
                    c_uptake = 1.0u"mol*mol^-1*hr^-1",
                ),
            ),  
            Params(
                assimilation_pars = ConstantNAssim(
                    n_uptake = 1.0u"mol*mol^-1*hr^-1",
                ),
            ),
        ),
        shared=shared,
        environment=env,
    )

    du2 = deepcopy(du0)
    u2 = deepcopy(u0)
    plant(du2, u2, nothing, 1u"hr")    

    # Carbon is balanced, with added C and N from assimilation
    @test du2[:VS] ≈ 1.0u"mol/hr" - du2[:CS] - maintenance_cost * 1u"mol" 
    @test du2[:VR] ≈ -du2[:CR] - maintenance_cost * 1u"mol"
    @test du2[:VS]/100 ≈ -du2[:NS] - maintenance_cost * 1u"mol"/100
    @test du2[:VR]/100 ≈ 1.0u"mol/hr" - du2[:NR] - maintenance_cost * 1u"mol"/100
end

@testset "balance with resorption is correct" begin
    plant3 = Plant(time = 0u"hr":1u"hr":8001u"hr", 
        params=(
            Params(
                assimilation_pars = ConstantCAssim(
                    c_uptake = 0.0u"mol*mol^-1*hr^-1",
                ),
            ),  
            Params(
                assimilation_pars = ConstantNAssim(
                    n_uptake = 0.0u"mol*mol^-1*hr^-1",
                ),
            ),
        ),
        shared=SharedParams(
            resorption_pars=LosslessResorption(),
            tempcorr_pars=nothing,
            su_pars = ParallelComplementarySU(),
            core_pars = DEBCore(
                j_E_mai = maintenance_cost,
                y_V_E = 1.0,
                y_E_C = 1.0,
                y_E_N = 100.0,
                n_N_V = 0.01,
                n_N_E = 0.01,
            ),
            catabolism_pars = CatabolismCNshared(
                k = 0.5u"d^-1",
            ),
        ),
        environment=env,
    )

    du3 = deepcopy(du0)
    u3 = deepcopy(u0)
    plant3(du3, u3, nothing, 1u"hr")    

    # Carbon is balanced, with losses to maintenance and returned resorption
    @test du3[:VS] == -du3[:CS] - maintenance_cost * 1u"mol"
    @test du3[:VR] == -du3[:CR] - maintenance_cost * 1u"mol"
    @test du3[:VS] ≈ -100du3[:NS] - maintenance_cost * 1u"mol"
    @test du3[:VR] ≈ -100du3[:NR] - maintenance_cost * 1u"mol"
end
