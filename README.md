# DEB

A Dynmaic Energy Budget modelling framework written in julia.

This is a generalised DEB model for plants, which can model organisms with N
structures and N reserves.

This model can be run in microclimates provided by the NicheMapr R
package - currently soil and air temperatures, soil moisture and insolation. 

It is also an in-progress attempt at using a functional programming style and
Julia's multiple-dispatch methods to abstract and generalise DEB theory and
maintain a short, maintainable codebase for multiple models.

The model.jl component could be used standalone but is best used in with MechanisticModels.jl
to provide easy setup, and interaction with solvers, optimisation and plotting tools.

Code adapted from original DEBtool plant model by Bas Kooijman:
https://github.com/add-my-pet/DEBtool_M
