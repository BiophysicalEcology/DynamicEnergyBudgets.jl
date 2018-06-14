module DynamicEnergyBudgets

using SimpleTraits
using AxisArrays
using StaticArrays
using DataFrames
using DataStructures
using Unitful
using Roots
using OrdinaryDiffEq
using Parameters
using Mixers
using MetaParameters
using NicheMap
using Photosynthesis

@metaparam label ""
@metaparam range [0.0, 1.0]

const TRANS = [:ass, :gro, :mai, :rep, :rej, :tra]
const TRANS1 = [:cat, :rej, :los]
const BI_XTOL = 1e-10u"d^-1"
const BI_MAXITER = 100

include("state.jl")
include("conversions.jl")
include("types.jl")
include("environment.jl")
include("assimilation.jl")
include("model.jl")
include("functions.jl")
include("apply.jl")
include("setup.jl")

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

export AbstractAssimilation,
       AbstractCarbonAssimilation,
       AbstractNitrogenAssimilation,
       NH4_NO3_Assimilation,
       N_Assimilation,
       C3Photosynthesis,
       KooijmanPhotosynthesis,
       KooijmanSLAPhotosynthesis,
       Kooijman_NH4_NO3_Assimilation,
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
       CarbonVars,
       NitrogenVars,
       Params,
       SharedParams,
       Vars,
       StateData,
       Organ,
       Organism,
       Scenario,
       OrganState,
       OrganismState,
       ScenarioState

end # module
