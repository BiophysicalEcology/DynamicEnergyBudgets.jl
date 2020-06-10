using LinearAlgebra,
      Unitful,
      LabelledArrays,
      DynamicEnergyBudgets,
      Test

using DynamicEnergyBudgets: translocate!, maintenence!, growth!, sum_flux!,
                            reserve_drain!, maturity!, metabolism!, catabolism!,
                            assimilation!, translocation!, shape, 
                            photosynthesis, define_organs, 
                            default, units

construct_organ(; params=Params(), shared=SharedParams(), vars=Vars()) =
    define_organs(Plant(params=(params,), vars=(vars,), shared=shared), 1)[1]


sh = Plant().shared
global n_ratios = [sh.n_N_P, sh.n_N_V, sh.n_N_M, sh.n_N_EC, sh.n_N_EN, sh.n_N_E]
# global n_ratios = [0.0, 0.14, 0.21, 0.01, 10.3, 0.24]

function factory()
    o = construct_organ()
    p = o.params
    u = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1e-2, 1e-2, 10.0]u"mol" 
    u = LVector{eltype(u),typeof(u),DynamicEnergyBudgets.STATE}(u)
    du = fill(0.0u"mol*hr^-1", 6)
    o, p, u, du
end

nothing
