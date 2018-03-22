module DynamicEnergyBudgets

using SimpleTraits
using AxisArrays
using StaticArrays
using DataFrames
using DataStructures
using MechanisticModels
using Photosynthesis
using NicheMap
using Plots

export DEBSettings, DEBStructure, DEBFlags, DEBFunctions
export @deb_settings
export integrate, ParamSpec, OrderedDict, DataFrame, Tspan 
export deb_widgets, deb_plottables

export StateV, StatePV, StateVE, StatePVE, StateVCN, StateVME, StatePVME, 
       StatePVCN, StatePVCNE, StatePVMCN, StatePVMCNE
export AbstractState, AbstractStateE, AbstractStateCN, AbstractStateCNE
export get_state1_names, init_state
export debplot


const FluxBase = AxisArray{Float64,3,Array{Float64,3},
                           Tuple{Axis{:state,Array{Symbol,1}},
                                 Axis{:transformations,Array{Symbol,1}},
                                 Axis{:time,UnitRange{Int64}}}}

const Flux = AxisArray{Float64,2,SubArray{Float64,2,Array{Float64,3},
                                          Tuple{Base.Slice{Base.OneTo{Int64}},
                                                Base.Slice{Base.OneTo{Int64}}, Int64},true},
                       Tuple{AxisArrays.Axis{:state,Array{Symbol,1}},
                             AxisArrays.Axis{:transformations,Array{Symbol,1}}}}

const StateBase = Vector{Vector{Float64}}

struct DEBFlags <: AbstractFlags
    exposed::FlagIndex
    time::FlagIndex
    temp::FlagIndex
end

struct DEBFunctions{F1,F2,F3,F4} <: AbstractFunctions
    area::F1
    assim::F2
    assim_sub::F3
    rate::F4
end

mutable struct DEBStructure{A,P,S,F<:DEBFunctions} <: AbstractStructure
    name::Symbol
    param_specs::ParamSpecs
    param_ids::Array{Int}
    init_params::P
    params::P
    flags::DEBFlags
    functions::F
    u::S
    A::Float64
    Jbase::FluxBase
    J1base::FluxBase
    J::Flux
    J1::Flux
    rates::Vector{Float64}
    assim_state::A
end

mutable struct DEBSettings{E,F,S<:Tuple{Vararg{<:DEBStructure}}} <: AbstractStructuredSettings
    u0::Vector{Float64}
    tspan::Tspan
    environment::E
    use_environment::Bool
    apply_environment!::F
    save_intermediate::Bool
    timestep_days::Float64
    structures::S
end

DEBSettings(u0, tspan, environment, use_environment, apply_environment!, 
            save_intermediate, timestep_days, state_type, structures) = begin
    DEBSettings(u0, tspan, environment, use_environment, apply_environment!, 
                save_intermediate, timestep_days, structures)
end

macro deb_settings(spec::Expr)
    block = MechanisticModels.main_block(spec)
    exp = settings_structured(:DEBSettings, block)
    return quote
        $exp
        apply(scale_time_dependent_params!, settings.structures, settings.timestep_days)
        apply(initialise_params!, settings.structures)
        settings
    end
end

include("conversions.jl")
include("constants.jl")
include("structures.jl")
include("state.jl")
include("assimilation_inputs.jl")
include("functions.jl")
include("model.jl")
include("environment.jl")
include("interface.jl")
include("plot.jl")

end # module
