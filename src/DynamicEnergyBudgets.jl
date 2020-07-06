"""
This package is a generalised DEB model. It was developed for plant modelling,
but can potentially be used to model any organisms and symbioses.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiplt models - potentially any organism.  Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M) plant model by Bas Kooijman.  """
module DynamicEnergyBudgets
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end DynamicEnergyBudgets

using ConstructionBase,
      DimensionalData,
      DocStringExtensions,
      FieldMetadata,
      FieldDefaults,
      FieldDocTables,
      Flatten,
      Mixers,
      OrdinaryDiffEq,
      Setfield,
      SimpleRoots,
      Requires,
      Unitful

using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Base: tail

import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable



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
       KooijmanNH4_NO3Assim, AbstractFvCBCAssim, BBPotentialCAssim, BallBerryCAssim, EmaxCAssim, ConstantCAssim

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

export AbstractScaling, Isomorph, V0morph, V1morph, V1V0morph, Plantmorph

export PassiveTranslocation, LosslessPassiveTranslocation, DissipativePassiveTranslocation

export ActiveTranslocation, LosslessActiveTranslocation, DissipativeActiveTranslocation

export AbstractCatabolism, Catabolism, CatabolismCN, CatabolismCNshared, CatabolismCNE

export AbstractAllometry, Allometry, SqrtAllometry, FixedAllometry

export AbstractParams, Params, SharedParams,
       Vars, PlottableVars, CarbonVars, NitrogenVars, Records, PlottableRecords

export AbstractOrgan, Organ, AbstractOrganism, Plant


const FIELDDOCTABLE = FieldDocTable((:Description, :Default, :Bounds),
                                    (description, default, bounds);
                                    truncation=(100,40,100))

const DEAD, ALIVE = false, true

# Auto docstrings
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    $(FIELDDOCTABLE)
    """

# Field metadata columns
@chain columns @udefault_kw @units @bounds @logscaled @description

const BI_XTOL = 1e-10
const BI_MAXITER = 100

include("traits.jl")
include("components/synthesizing_units.jl")
include("components/temperature_correction.jl")
include("components/scaling.jl")
include("components/allometry.jl")
include("components/resorption.jl")
include("components/assimilation.jl")
include("components/germination.jl")
include("components/maturity.jl")
include("components/production.jl")
include("components/translocation.jl")
include("components/core.jl")
include("components/catabolism.jl")
include("organism.jl")
include("components/environment.jl")
include("functions.jl")
include("setup.jl")
include("model.jl")

function __init__()
    @require Photosynthesis="5e10c064-2706-53a3-a67d-d473e313a663" include("components/fvcb.jl")
    @require Microclimate="e0eb800d-4a9f-54ae-b0f8-217228d9d7c3" include("components/microclimate.jl")
end

end # module
