"""
This is a generalised DEB model. It was developed for plant modelling, but can potentially
model any organisms and symbioses.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiplt models - potentially any organism.  Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M) plant model by Bas Kooijman.  """ 
module DynamicEnergyBudgets

using Unitful,
      OrdinaryDiffEq,
      DocStringExtensions,
      MacroTools,
      Distributions,
      LabelledArrays,
      SimpleRoots,
      Mixers,
      FieldMetadata,
      FieldDefaults,
      Photosynthesis,
      Microclimate,
      Apply,
      Flatten

using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Base: tail

import FieldDefaults: get_default
import FieldMetadata: @prior, @default, @redefault, @description, @units, @limits, @logscaled, @flattenable, @plottable, @selectable,
                      prior, default, description, units, limits, logscaled, flattenable, plottable, selectable
import Photosynthesis: potential_dependence                      


export tempcorr,
       rate_bracket,
       rate_formula,
       catabolic_fluxes,
       half_saturation,
       stoich_merge,
       synthesizing_unit

export integrate,
       get_state1_names,
       init_state,
       describe,
       paramrange,
       runmodel!,
       debmodel!,
       apply_environment!

export AbstractAssimilation,
       AbstractNAssim, NH4_NO3Assim, NAssim, ConstantNAssim,
       AbstractCAssim, C3Photosynthesis,
       KooijmanPhotosynthesis, KooijmanSLAPhotosynthesis, KooijmanWaterPotentialPhotosynthesis, 
       KooijmanNH4_NO3Assim, AbstractFvCBCAssim, BallBerryCAssim, EmaxCAssim, ConstantCAssim

export AbstractSynthesizingUnit, ParallelComplementarySU, MinimumRuleSU, KfamilySU

export AbstractRate, SimpleRate, FZeroRate

export AbstractResorption, LosslessResorption, StructuralLossResorption, DissipativeResorption

export AbstractTemperatureCorrection, TempCorr, TempCorrLower, TempCorrLowerUpper, ParentTardieu

export AbstractScaling, KooijmanArea

export AbstractMaturity, Maturity

export AbstractMaintenance, Maintenance

export AbstractDEBCore, DEBCore

export AbstractGermination, ThresholdGermination

export AbstractProduction, Production

export AbstractShape, Isomorph, V0morph, V1morph, V1V0morph, Plantmorph

export AbstractRejection, LosslessRejection, DissipativeRejection

export AbstractCatabolism, Catabolism, CatabolismCN, CatabolismCNshared, CatabolismCNE

export AbstractTranslocation, LosslessMultipleTranslocation, DissipativeMultipleTranslocation,
       LosslessTranslocation, DissipativeTranslocation

export AbstractAllometry, Allometry, SqrtAllometry, FixedAllometry

export AbstractParams, Params, SharedParams,
       Vars, CarbonVars, NitrogenVars, Records,
       AbstractOrgan, Organ,
       AbstractOrganism, Plant


# Auto docstrings
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

# Field metadata columns
@chain columns @description @logscaled @limits @prior @units @udefault_kw

const STATELEN = 6

const STATE = (:P, :V, :M, :C, :N, :E)
const STATE1 = (:EE, :CN, :C, :N, :E)
const TRANS = (:ass, :gro, :mai, :mat, :rej, :tra, :fbk)
const TRANS1 = (:ctb, :rej, :los)
const BI_XTOL = 1e-10
const BI_MAXITER = 100

include("traits.jl")
include("components/synthesizing_units.jl")
include("components/temperature_correction.jl")
include("components/rate.jl")
include("components/shape.jl")
include("components/allometry.jl")
include("components/resorption.jl")
include("components/assimilation.jl")
include("components/germination.jl")
include("components/maturity.jl")
include("components/production.jl")
include("components/translocation.jl")
include("components/catabolism.jl")
include("components/maintenance.jl")
include("types.jl")
include("environment.jl")
include("getters.jl")
include("functions.jl")
include("setup.jl")
include("model.jl")
include("aliases.jl")

end # module
