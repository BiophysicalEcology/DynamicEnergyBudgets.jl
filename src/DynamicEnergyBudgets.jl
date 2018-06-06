module DynamicEnergyBudgets

using SimpleTraits
using AxisArrays
using StaticArrays
using DataFrames
using DataStructures
using Unitful
using SimpleRoots
using OrdinaryDiffEq
using Parameters
using Mixers
using MetaParameters
using NicheMap
using Photosynthesis

@metaparam label ""
@metaparam range [0.0, 1.0]

include("state.jl")
include("constants.jl")
include("conversions.jl")
include("types.jl")
include("environment.jl")
include("assimilation_inputs.jl")
include("model.jl")
include("functions.jl")
include("apply.jl")

export tempcorr,
       rate_bracket,
       rate_formula,
       catabolic_fluxes,
       half_saturation,
       stoich_merge,
       synthesizing_unit

export integrate,
       Tspan,
       get_state1_names,
       init_state,
       describe,
       paramrange,
       runmodel!,
       debmodel!,
       apply_environment!

export AbstractState, 
       AbstractStateE, 
       AbstractStateCN, 
       AbstractStateCNE,
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
       NitrogenVars,
       Params,
       Vars,
       StateData,
       Organ,
       Organism,
       Scenario,
       OrganState,
       OrganismState,
       ScenarioState


end # module
