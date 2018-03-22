using Revise

using Base.Test
using StaticArrays
using DataStructures
using DataFrames
using MechanisticModels
using DynamicEnergyBudgets
using MechanisticModels.cascade_update_params!


uptv!(tv) = begin
    tp.param_three = 99.0
end

tspan = (1.0, 10.0)
MechanisticModels.eval_paramstype([:param_one, :param_two, :param_three])
ps1 = MechanisticModels.ParamSpecs(:param_one => ParamSpec(1, 1.0), :param_two => ParamSpec(2, 2.0, :time))
ps2 = MechanisticModels.ParamSpecs(:param_two => ParamSpec(2, 99.0), :param_three => ParamSpec(3, 3.0))

function build_test_structure(param_specs)
    v = MechanisticModels._Parameters(0.0, 0.0, 0.0)
    fl = MechanisticModels.get_flags(DEBStructure, param_specs)
    f = DEBFunctions(x -> x, x -> 2, x -> 3x, x -> 4x)
    settings = OrderedDict(:u0 => [0.0, 0.0, 0.0, 0.0], :state_type => StateVE, :tspan => tspan)
    DEBStructure(:name, param_specs, v, fl, f, settings)
end


@testset "update structure" begin
    structures = build_test_structure.((ps1, ps2)) 
    s1 = structures[1]
    s2 = structures[2]
    u = zeros(4)
    du = zeros(4)
    num_structures = 2
    states = 2
    trans = length(DynamicEnergyBudgets.TRANS)
    timestep_days = 1.0/24.0
    offset = 0

    @testset "sum_flux! sums all J arrays correctly" begin
        fill!(structures[1].J, 1.0)
        DynamicEnergyBudgets.sum_flux!(du, structures, 2, offset, states, trans)
        @test du == [6.0, 6.0, 0.0, 0.0]                                       
        fill!(structures[2].J, 1.0)                                        
        DynamicEnergyBudgets.sum_flux!(du, structures, 2, offset, states, trans)
        @test du == [6.0, 6.0, 6.0, 6.0]                                       
    end

    @testset "split_state! copies state into structures" begin
        structures = build_test_structure.((ps1, ps2)) 
        @test structures[1].u == StateVE(0.0, 0.0)
        @test structures[2].u == StateVE(0.0, 0.0)

        u = [1.0, 2.0, 3.0, 4.0]
        DynamicEnergyBudgets.split_state!(structures, 0, u)
        @test structures[1].u == StateVE(1.0, 2.0)
        @test structures[2].u == StateVE(3.0, 4.0)
    end

    @testset "set_current_flux! assigns J to correct view of Jbase" begin
        structures = build_test_structure.((ps1,  ps2)) 
        s1 = structures[1]
        s2 = structures[2]
        DynamicEnergyBudgets.set_current_flux!(structures, 0)
        @test s1.Jbase[1, 1, 1] == 0.0
        s1.J[1, 1, 1] =  100.0
        @test s1.Jbase[1, 1, 1] == 100.0
        @test s1.Jbase[1, 1, 2] == 0.0
    end

    @testset "cascade update and initialise params" begin
        structures = build_test_structure.((ps1, ps2)) 
        @test structures[1].params.param_one == 0.0
        @test structures[1].param_specs[:param_one].value == 1.0
        @test_throws KeyError structures[2].param_specs[:param_one].value
        cascade_update_params!(structures)
        DynamicEnergyBudgets.initialise_params!(structures)
        # param_one cascades down to structure 2
        @test structures[1].params.param_one == 1.0
        @test structures[2].params.param_one == 1.0
        # Both update param_two
        @test structures[1].params.param_two == 2.0
        @test structures[2].params.param_two == 99.0
        # param_three does not cascade back to structure 1
        @test structures[1].params.param_three == 0.0
        @test structures[2].params.param_three == 3.0
    end

    @testset "scale_time_dependent_params! scales :time flagged params" begin
        structures = build_test_structure.((ps1, ps2)) 
        cascade_update_params!(structures)
        @test structures[2].init_params.param_three == 3.0
        DynamicEnergyBudgets.scale_time_dependent_params!(structures, timestep_days)
        # params with :time flag are scaled
        @test structures[1].init_params.param_two == timestep_days * 2.0 
        # Others remain unaffected
        @test structures[2].init_params.param_two == 99.0 
        @test structures[1].init_params.param_one == 1.0 
    end
end

