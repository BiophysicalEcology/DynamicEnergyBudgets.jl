# DynamicEnergyBudgets

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:module]
```

## Index

### Types

```@index
Modules = [DynamicEnergyBudgets]
Order   = [:type]
```

### Functions

```@index
Modules = [DynamicEnergyBudgets]
Order   = [:function]
```

## Parameters and Types

DynamicEnergyBudgets uses the julia type system to organise 
both the parameters and the logic of the model. This means the model is
completely modular and customisable.

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:type]
Pages   = ["types.jl"]
```

## Model

Core functions of the model.

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:function]
Pages   = ["model.jl"]
```

## Environment

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:function]
Pages   = ["environment.jl"]
```

## Assimilation

Carbon and nitrogen assimilation functions. These are the main points of
difference between organs.

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:function]
Pages   = ["assimilation.jl"]
```

## Functions

Low-level DEB theory functions

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:function]
Pages   = ["functions.jl"]
```

## Setup and utilities

```@autodocs
Modules = [DynamicEnergyBudgets]
Order   = [:function]
Pages   = ["apply.jl"]
Pages   = ["setup.jl"]
```
