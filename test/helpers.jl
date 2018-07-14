using DynamicEnergyBudgets: build_J, build_J1

construct_organ(; params=Params(), shared=SharedParams(), vars=Vars(),
      J=build_J(1.0u"mol/hr", typeof(1.0u"mol/hr")),
      J1=build_J1(1.0u"mol/hr", typeof(1.0u"mol/hr"))) = begin
    Organ(params, shared, vars, J, J1)
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
    u1 = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    u2 = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    u2 .*= 2.7
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    o1.vars.scale = scaling(o1.params.scaling, u1[V])
    o2.vars.scale = scaling(o2.params.scaling, u2[V])
    f = N_Assimilation();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end

function cfactory()
    o1 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8), 
               vars=Vars(assimilation=CarbonVars()));;
    p1 = o1.params;
    o2 = construct_organ(params=Params(y_EC_ECT = 0.8, y_EN_ENT = 0.8));
    p2 = o2.params;
    u1 = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    u2 = [9.0u"mol",8.0u"mol",7.0u"mol",6.0u"mol",5.0u"mol",4.0u"mol"]
    u2 .*= 1.9
    du1 = fill(0.0u"mol*hr^-1", 6)
    du2 = fill(0.0u"mol*hr^-1", 6)
    o1.vars.scale = scaling(o1.params.scaling, u1[V])
    o2.vars.scale = scaling(o2.params.scaling, u2[V])
    f = KooijmanSLAPhotosynthesis();

    o1, p1, u1, du1, o2, p2, u2, du2, f
end
