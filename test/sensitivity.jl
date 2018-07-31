using Revise
using DynamicEnergyBudgets
using CompositeFieldVectors
using AxisArrays
using ForwardDiff
using DiffEqSensitivity
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

du = [0.0 for i in 1:12]
u = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1e-2, 1e-2, 10.0]
organism = DynamicEnergyBudgets.UntypedPlant();
p = [flatten(organism)...];
n = [flatten_fieldnames(organism)...];
np = AxisArray(p, Axis{:parameters}(n));
jac = ForwardDiff.jacobian(p1->organism(du, u, p1, 1), np)
@test sum(jac) != 0.0

organism.records[1].J
