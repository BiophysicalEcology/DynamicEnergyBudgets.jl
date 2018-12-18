using Revise,
      LinearAlgebra,
      Unitful,
      LabelledArrays,
      DynamicEnergyBudgets,
      Test

using DynamicEnergyBudgets: reuse_rejected!, translocate!, maintenence!, growth!, sum_flux!,
                            reserve_drain!, reserve_loss!, maturity!, metabolism!, catabolism!,
                            assimilation!, translocation!, scaling, set_scaling!,
                            uptake_nitrogen, photosynthesis, STATELEN,
                            define_organs, default, units, set_var!, rate, unpack, tempcorrection, θE

construct_organ(; params=Params(), shared=SharedParams(), vars=Vars()) = begin
    define_organs(Organism(params=(params,), vars=(vars,), shared=shared), 1)[1]
end


sh = Organism().shared
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
