# DynamicEnergyBudgets

[![Build Status](https://travis-ci.org/rafaqz/DynamicEnergyBudgets.jl.svg?branch=master)](https://travis-ci.org/rafaqz/DynamicEnergyBudgets.jl)
[![Coverage Status](https://coveralls.io/repos/rafaqz/DynamicEnergyBudgets.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/rafaqz/DynamicEnergyBudgets.jl?branch=master)
[![codecov.io](http://codecov.io/github/rafaqz/DynamicEnergyBudgets.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Microclimate.jl?branch=master)

A Dynamic Energy Budget modelling framework written in Julia.

This is a generalised DEB model that for plant modelling, but can
used to used model any kind of organism.

This models can also be run in microclimates provided by the NicheMapR R package
using [Microclimate.jl](https://github.com/rafaqz/Microclimate.jl), 
and can use wide a range of photosynthesis and stomatal conductance formulations 
from [Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).


Code is largely adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M)
plant model by Bas Kooijman.
