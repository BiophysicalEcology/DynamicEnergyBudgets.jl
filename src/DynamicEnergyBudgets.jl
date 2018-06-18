"""
This is a generalised DEB model. It was developed for plant modelling, but can potentially 
model any organisms and symbioses.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to
abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiple models - potentially any organism.

Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M)
plant model by Bas Kooijman.
"""
module DynamicEnergyBudgets

using SimpleTraits,
      AxisArrays,
      StaticArrays,
      DataFrames,
      Unitful,
      Roots,
      OrdinaryDiffEq,
      Parameters,
      Mixers,
      MetaParameters,
      Microclimate,
      Photosynthesis,
      DocStringExtensions

@metaparam label ""
@metaparam range [0.0, 1.0]

@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

Base.muladd(a::Quantity, b::Quantity, c::Quantity) = a * b + c

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
