# DynamicEnergyBudgets

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:module]
```

## Model

Basic types functions of a DEB model.

```@docs
AbstractOrgan
Organ
AbstractVars
Vars
PlottableVars
AbstractParams
Params
AbstractSharedParams
SharedParams
AbstractOrganism
Plant
DynamicEnergyBudgets.debmodel!
DynamicEnergyBudgets.metabolism!
```


## Components

### Core Parameters

```@docs
DEBCore
DynamicEnergyBudgets.growth!
DynamicEnergyBudgets.maintenance!
```

### Allometry

```@docs
AbstractAllometry
Allometry
SqrtAllometry
FixedAllometry
```

### Assimilation

```@docs
AbstractAssim
AbstractCAssim
ConstantCAssim
KooijmanSLAPhotosynthesis
KooijmanWaterPotentialPhotosynthesis
BallBerryCAssim
BallBerryPotentialCAssim
CarbonVars
AbstractNAssim
ConstantNAssim
NAssim
KooijmanNH4_NO3Assim
NitrogenVars
DynamicEnergyBudgets.assimilation!
DynamicEnergyBudgets.photosynthesis
DynamicEnergyBudgets.nitrogen_uptake
```

### Catabolism

```@docs
CatabolismCN
CatabolismCNE
CatabolismCNshared
DynamicEnergyBudgets.catabolism!
```

### Environment

```@docs
ManualTemperature
DynamicEnergyBudgets.apply_environment!
```

### Germination

```@docs
AbstractGermination
ThresholdGermination
DynamicEnergyBudgets.isgerminated
```

### Maturity

```@docs
AbstractMaturity
Maturity
DynamicEnergyBudgets.maturity!
```

### Production

```@docs
Production
```

### Rate

```@docs
DynamicEnergyBudgets.calc_rate
DynamicEnergyBudgets.rate_formula
```

### Resorption

```@docs
AbstractResorption
LosslessResorption
StructuralLossResorption
DissipativeResorption
DynamicEnergyBudgets.resorption
```

### Scaling

```@docs
AbstractScaling
Isomorph
V0morph
V1morph
V1V0morph
Plantmorph
DynamicEnergyBudgets.scaling_correction
```

### Synthesizing Units

```@docs
AbstractSynthesizingUnit
ParallelComplementarySU
MinimumRuleSU
KfamilySU
DynamicEnergyBudgets.synthesizing_unit
DynamicEnergyBudgets.stoich_merge
```

### Temperature Correction

```@docs
AbstractTemperatureCorrection
TempCorr
TempCorrLower
TempCorrLowerUpper
ParentTardieu
DynamicEnergyBudgets.tempcorr
```

### Translocation

```@docs
DissipativePassiveTranslocation
LosslessPassiveTranslocation
LosslessActiveTranslocation
DissipativeActiveTranslocation
DynamicEnergyBudgets.active_translocation!
DynamicEnergyBudgets.passive_translocation!
DynamicEnergyBudgets.translocation!
```

## Other functions

Low-level model functions

```@autodocs
Modules = [DynamicEnergyBudgets]
Pages   = ["functions.jl", "setup.jl"]
```
