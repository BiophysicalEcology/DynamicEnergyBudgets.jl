using Revise
using DynamicEnergyBudgets
using DynamicEnergyBudgets: reuse_rejected!, dissipation!, translocate!, product!, 
                            maintenence!, growth!, sum_flux!, reserve_drain!, reserve_loss!,
                            maturity!, metabolism!, catabolism!, assimilation!, translocation!,
                            scaling, P, V, M, C, N, E, EE, CN, STATELEN, ass, gro, mai, rep, rej, tra, cat, rej, los,

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("helpers.jl")
include("setup.jl")
include("math.jl")
include("assimilation.jl")
include("environment.jl")
include("balance.jl")
include("diffeq.jl")
