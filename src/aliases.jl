FLUX = 1.0u"mol/hr"

ShootParams(;kwargs...) = Params(name=:shoot, assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
RootParams(;kwargs...) = Params(name=:root, assimilation=N_Assimilation(), kwargs...)
ShootVars(;kwargs...) = Vars(assimilation=CarbonVars(), kwargs...)
RootVars(;kwargs...) = Vars(assimilation=NitrogenVars(), kwargs...)
Records(params; vars=(ShootVars(), RootVars()), t=1:1:1, val=FLUX) = begin
    (Records(params[1], vars[1], t, val, typeof(val)), 
     Records(params[2], vars[2], t, val, typeof(val)))
end
Plant(;kwargs...) = Organism(;params=(ShootParams(), RootParams()), vars=(ShootVars(), RootVars()), kwargs...)

UntypedShootVars(;kwargs...) = 
    default_kw(Vars{Any,Any,Any,Any,Any}; assimilation=default_kw(CarbonVars{Any,Any}), kwargs...)
UntypedRootVars(;kwargs...) = 
    default_kw(Vars{Any,Any,Any,Any,Any}; assimilation=default_kw(NitrogenVars{Any,Any}), kwargs...)
UntypedRecords(params; vars=(UntypedShootVars(), UntypedRootVars()), t=1:1:1, val=FLUX, typ=ANY) = begin
    (Records(params[1], vars[1], t, FLUX, typ), 
     Records(params[2], vars[2], t, FLUX, typ))
end
"Outer construtor for defaults"
UntypedPlant(; params = (ShootParams(), RootParams()),
           vars = (UntypedShootVars(), UntypedRootVars()),
           shared = SharedParams(),
           environment = nothing,                       
           time = 0u"hr":1u"hr":1000u"hr") = begin
    records = UntypedRecords(params; vars=vars)
    Organism(params, shared, records, environment)
end

FvCBShootParams(;kwargs...) = Params(name=:shoot, assimilation=FvCBPhotosynthesis(), kwargs...)
FvCBRootParams(;kwargs...) = Params(name=:root, assimilation=N_Assimilation(), kwargs...)
FvCBShootVars(;kwargs...) = Vars(assimilation=Photosynthesis.PhotoVars(), kwargs...)
FvCBRootVars(;kwargs...) = Vars(;assimilation=NitrogenVars(), kwargs...)
FvCBPlant(;kwargs...) = Organism(;params=(FvCBShootParams(), FvCBRootParams()),
                                  vars=(FvCBShootVars(), FvCBRootVars()), kwargs...)
FvCBNoRoots(;kwargs...) = Organism(;params=(FvCBShoot()), kwargs...)

ConstantShoot(;kwargs...) = Organ(name=:shoot, params=Params(assimilation=ConstantCarbonAssimilation()), kwargs...)
ConstantRoot(;kwargs...) = Organ(name=:root, params=Params(assimilation=ConstantNitrogenAssimilation()), kwargs...)
ConstantPlant(;kwargs...) = Organism(; params=(ConstantShootParams(), ConstantRootParams()), kwargs...)

LeafParams(;kwargs...) = Params(name=:leaf, translocation=ShootTrans(dest=:stem), 
                                assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
StemParams(;kwargs...) = Params(name=:stem, translocation=RootTrans(dest=(:leaf, :root), prop=0.5), 
                                assimilation=N_Assimilation(), kwargs...)
RootParams3(;kwargs...) = Params(name=:root, translocation=RootTrans(dest=:stem), 
                                 assimilation=N_Assimilation(), kwargs...)
