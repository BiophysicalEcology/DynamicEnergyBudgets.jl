"""
This is a generalised DEB model. It was developed for plant modelling, but can potentially 
model any organisms and symbioses.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiplt models - potentially any organism.  Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M) plant model by Bas Kooijman.  """ module DynamicEnergyBudgets 
using Unitful,
      OrdinaryDiffEq,
      ForwardDiff,
      DocStringExtensions,
      Distributions,
      SimpleRoots,
      Mixers,
      Tags,
      Defaults,
      Microclimate,
      Photosynthesis,
      Flatten

import Flatten: flattenable
import Defaults: get_default
import Tags: @prior, @default, @description, @units, @limits, prior, default, description, units, limits

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
       Params,
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

const P = 1
const V = 2
const M = 3
const C = 4
const N = 5
const E = 6

const EE = 1
const CN = 2

const STATELEN = 6

const ass = 1
const gro = 2 
const rej = 3 
const mai = 4 
const mat = 5 
const tra = 6
const fbk = 7

const cat = 1
const los = 2
const rej = 3

const STATE = [:P, :V, :M, :C, :N, :E]
const STATE1 = [:EE, :CN, :_, :C, :N, :E]
const TRANS = [:ass, :gro, :mai, :mat, :rej, :tra, :fbk]
const TRANS1 = [:cat, :rej, :los]
const BI_XTOL = 1e-10
const BI_MAXITER = 100

# include("state.jl")
include("types.jl")
include("aliases.jl")
include("environment.jl")
include("assimilation.jl")
include("model.jl")
include("functions.jl")
include("apply.jl")
include("setup.jl")
include("sensitivity.jl")

end # module
