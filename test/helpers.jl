using Revise, 
      Unitful, 
      DynamicEnergyBudgets,
      Test

using DynamicEnergyBudgets: reuse_rejected!, dissipation!, translocate!, product!, 
                            maintenence!, growth!, sum_flux!, reserve_drain!, reserve_loss!,
                            maturity!, metabolism!, catabolism!, assimilation!, translocation!,
                            scaling, P, V, M, C, N, E, EE, CN, STATELEN, ass, gro, mai, mat, rej, tra, cat, rej, los, 
                            define_organs, default, units, set_var!

construct_organ(; params=Params(), shared=SharedParams(), vars=Vars()) = begin
    define_organs(Organism(params=(params,), vars=(vars,), shared=shared), 1)[1]
end

function factory()
    o = construct_organ()
    p = o.params
    u = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    du = fill(0.0u"mol*hr^-1", 6)
    o, p, u, du
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
