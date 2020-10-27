# DynamicEnergyBudgets

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/DynamicEnergyBudgets.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/DynamicEnergyBudgets.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/DynamicEnergyBudgets.jl.svg?branch=master)](https://travis-ci.com/rafaqz/DynamicEnergyBudgets.jl)
[![codecov.io](http://codecov.io/github/rafaqz/DynamicEnergyBudgets.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/DynamicEnergyBudgets.jl?branch=master)

A Dynamic Energy Budget modelling framework written in Julia.

This is a generalised DEB model that for plant modelling, but can
used to used model any kind of organism.

This models can also be run in microclimates provided by the NicheMapR R package
using [Microclimate.jl](https://github.com/rafaqz/Microclimate.jl), 
and can use wide a range of photosynthesis and stomatal conductance formulations 
from [Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

See scripts at https://github.com/rafaqz/DEBplant for a live user interface and plotting examples. As this package and Photosynthesis.jl are officially registered, so build a project using DEBplant as the package versions are locked in the Manifest.toml.

If you need this package to be registered, make an issue and request it! I am not currently working on DEB models and have many other packages competing for my time. If there is a project that wishes to use it, I will gladly help get things working for you. 


Code is largely adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M)
plant model by Bas Kooijman.


![Plant model](https://raw.githubusercontent.com/rafaqz/DynamicEnergyBudgets.jl/assets/deb_plant.png)
