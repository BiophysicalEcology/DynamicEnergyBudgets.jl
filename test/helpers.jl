using Revise,
      LinearAlgebra,
      Unitful,
      DynamicEnergyBudgets,
      Test

using DynamicEnergyBudgets: reuse_rejected!, translocate!, maintenence!, growth!, sum_flux!,
                            reserve_drain!, reserve_loss!, maturity!, metabolism!, catabolism!,
                            assimilation!, translocation!, set_scaling!, P, V, M, C, N, E, EE, CN,
                            STATELEN, ass, gro, mai, mat, rej, tra, ctb, rej, los,
                            define_organs, default, units, set_var!, rate, unpack, tempcorrection, θE

construct_organ(; params=Params(), shared=SharedParams(), vars=Vars()) = begin
    define_organs(Organism(params=(params,), vars=(vars,), shared=shared), 1)[1]
end

function factory()
    o = construct_organ()
    p = o.params
    u = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1e-2, 1e-2, 10.0]u"mol"
    du = fill(0.0u"mol*hr^-1", 6)
    sh = o.shared
    n_ratios = [sh.n_N_P, sh.n_N_V, sh.n_N_M, sh.n_N_EC, sh.n_N_EN, sh.n_N_E]
    o, p, u, du, n_ratios
end

function nfactory()
    o1 = construct_organ(params=Params(y_EC_ECT=0.8, y_EN_ENT=0.8),
               vars=Vars(assimilation=NitrogenVars()));;
    p1 = o1.params
    o2 = construct_organ(params=Params(y_EC_ECT=0.8, y_EN_ENT=0.8))
    p2 = o2.params;
    u1 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 .*= 2.7
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    set_var!(o1.vars, :scale, scaling(o1.params.scaling, u1[V]))
    set_var!(o2.vars, :scale, scaling(o2.params.scaling, u2[V]))
    f = NAssim();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

function cfactory()
    o1 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8),
               vars=Vars(assimilation=CarbonVars()));;
    p1 = o1.params;
    o2 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8));
    p2 = o2.params;
    u1 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 = [9.0,8.0,7.0,6.0,5.0,4.0]u"mol"
    u2 .*= 1.9
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    set_var!(o1.vars, :scale, scaling(o1.params.scaling, u1[V]))
    set_var!(o2.vars, :scale, scaling(o2.params.scaling, u2[V]))
    f = KooijmanSLAPhotosynthesis();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

nothing
