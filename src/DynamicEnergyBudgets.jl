module DynamicEnergyBudgets

using SimpleTraits
using AxisArrays
using StaticArrays
using DataFrames
using DataStructures
using Unitful
using SimpleRoots
using MultiScaleArrays
using DifferentialEquations
using Parameters

using UnitfulMoles
# using BiophysicalModels
using DynamicEnergyBudgetsBase
using PlantPhysiology
using NicheMap


macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

include("state.jl")
include("constants.jl")
include("conversions.jl")
include("types.jl")
include("environment.jl")
# include("structures.jl")
include("assimilation_inputs.jl")
include("model.jl")


export integrate, OrderedDict, DataFrame, Tspan 

export get_state1_names, init_state
export area_mass_kooijman, shoot_assimilation!, maestra,
       photosynthesis_sla, photosynthesis_kooijman, photomoduleC3,
       root_assimilation!, nitrogen_uptake, nitrogen_uptake_sla
export deb_model, find_rate
export load_nichemap, apply_nichemap!, apply_deb_environment!, apply_energy_balance!, apply_maestra!

export AbstractState, 
       AbstractStateE, 
       AbstractStateCN, 
       AbstractStateCNE
       StateV,
       StateE,
       StatePV,
       StateVE,
       StateCN,
       StatePVE,
       StateVCN,
       StateVME,
       StateCNE,
       StatePVME,
       StatePVCN,
       StatePVMCN,
       StatePVCNE,
       StatePVMCNE

export Assimilation,
       CarbonAssimilation,
       NitrogenAssimilation,
       NH4_NO3_Assimilation,
       KooijmanPhotosynthesis,
       C3Photosynthesis,
       KooijmanSLAPhotosynthesis,
       Kooijman_NH4_NO3_Assimilation,
       KooijmanSLA_NH4_NO3_Assimilation,
       AbstractStateFeedback,
       Autophagy,
       AbstractTempCorr,
       TempCorr,
       TempCorrLower,
       TempCorrLowerUpper,
       AbstractScaling,
       KooijmanArea,
       Structure,
       Products,
       Maturity,
       CarbonReserve,
       NitrogenReserve,
       GeneralReserve,
       Params,
       StateData,
       Organ,
       Organism,
       Scenario

end # module
