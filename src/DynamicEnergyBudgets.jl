"""
This is a generalised DEB model. It was developed for plant modelling, but can potentially 
model any organisms and symbioses.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiplt models - potentially any organism.  Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M) plant model by Bas Kooijman.  """ module DynamicEnergyBudgets 
using LinearAlgebra,
      Unitful,
      OrdinaryDiffEq,
      ForwardDiff,
      DocStringExtensions,
      Distributions,
      SimpleRoots,
      Mixers,
      LabelledArrays,
      FieldMetadata,
      Defaults,
      # Microclimate,
      Photosynthesis,
      UnitlessFlatten

using Base: tail
import UnitlessFlatten: flattenable
import Defaults: get_default
import FieldMetadata: @prior, @default, @description, @units, @limits, prior, default, description, units, limits

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

export AbstractAssim,
       AbstractCAssim,
       AbstractNAssim,
       NH4_NO3Assim,
       NAssim,
       C3Photosynthesis,
       KooijmanPhotosynthesis,
       KooijmanSLAPhotosynthesis,
       KooijmanNH4_NO3Assim,
       ConstantCAssim,
       ConstantNAssim,
       AbstractStateFeedback,
       Autophagy,
       AbstractTempCorr,
       TempCorr,
       TempCorrLower,
       TempCorrLowerUpper,
       AbstractScaling,
       KooijmanArea,
       Maturity,
       CarbonVars,
       NitrogenVars,
       ParamsCNE,
       ParamsCN,
       SharedParams,
       Vars,
       Organ,
       Organism,
       Records

@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

const STATELEN = 6

const STATE = (:P, :V, :M, :C, :N, :E)
const STATE1 = (:EE, :CN, :C, :N, :E)
const TRANS = (:ass, :gro, :mai, :mat, :rej, :tra, :fbk)
const TRANS1 = (:ctb, :rej, :los)
const BI_XTOL = 1e-10
const BI_MAXITER = 100

include("types.jl")
include("traits.jl")
include("aliases.jl")
# include("environment.jl")
include("assimilation.jl")
include("model.jl")
include("functions.jl")
include("apply.jl")
include("setup.jl")
include("sensitivity.jl")

end # module
