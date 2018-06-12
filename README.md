# DynamicEnergyBudgets

A Dynmaic Energy Budget modelling framework written in Julia.

This is a generalised DEB model that for plant modelling, but can potentially 
model organisms and symbioses with N organs and N reserves.

This model can also be run in microclimates provided by the NicheMapr R package, and
can use wide a range of photosynthesis and stomatal conductance formulations from
[Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

It is also an in-progress attempt at using Julia's multiple-dispatch methods to
abstract and generalise DEB theory and maintain a short, maintainable codebase
for multiple models - potentially any organism.

Code is adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M)
plant model by Bas Kooijman.
