FLUX = 1.0u"mol/hr"

ShootParams(;kwargs...) = Params(;name=:shoot, assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
RootParams(;kwargs...) = Params(;name=:root, assimilation=NAssim(), kwargs...)
ShootVars(;kwargs...) = Vars(;assimilation=CarbonVars(), kwargs...)
RootVars(;kwargs...) = Vars(;assimilation=NitrogenVars(), kwargs...)
Records(params; vars=(ShootVars(), RootVars()), t=1:1:1, val=FLUX) = begin
    (Records(params[1], vars[1], t, val, typeof(val)), 
     Records(params[2], vars[2], t, val, typeof(val)))
end
Plant(;kwargs...) = Organism(;params=(ShootParams(), RootParams()), vars=(ShootVars(), RootVars()), kwargs...)


NoMaturityShootParams(;kwargs...) = ShootParams(; maturity=nothing, kwargs...)
NoMaturityRootParams(;kwargs...) = RootParams(;maturity=nothing, kwargs...)
NoMaturityPlant(;kwargs...) = Organism(;params=(NoMaturityShootParams(), NoMaturityRootParams()), vars=(ShootVars(), RootVars()), kwargs...)


FvCBShootParams(;kwargs...) = Params(;name=:shoot, assimilation=FvCBPhotosynthesis(), kwargs...)
FvCBLeafParams(;kwargs...) = Params(;name=:shoot, assimilation=FvCBPhotosynthesis(),
                                    translocation=Translocation(destnames=:stem), kwargs...)
FvCBStemParams(;kwargs...) = Params(;name=:shoot, assimilation=nothing, 
                                    translocation=Translocation((:leaf,:root), 0.5), kwargs...)
FvCBRootParams(;kwargs...) = Params(;name=:root, translocation=Translocation(destnames=:stem), 
                                    assimilation=NAssim(), kwargs...)
FvCBShootVars(;kwargs...) = Vars(;assimilation=Photosynthesis.PhotoVars(), kwargs...)
FvCBStemVars(;kwargs...) = Vars(;assimilation=nothing, kwargs...)
FvCBLeafVars(;kwargs...) = Vars(;assimilation=Photosynthesis.PhotoVars(), kwargs...)
FvCBRootVars(;kwargs...) = Vars(;assimilation=NitrogenVars(), kwargs...)
FvCBNoRoots(;kwargs...) = Organism(;params=(FvCBShoot()), kwargs...)
FvCBPlant(;kwargs...) = Organism(;params=(FvCBShootParams(), FvCBRootParams()),
                                  vars=(FvCBShootVars(), FvCBRootVars()), kwargs...)
FvCBPlant3(;kwargs...) = Organism(;params=(FvCBLeafParams(), FvCBStemParams(), FvCBRootParams()), 
                              vars=(FvCBLeafVars(), FvCBStemVars(), FvCBRootVars()), kwargs...)


ConstantShootParams(;kwargs...) = Params(assimilation=ConstantCAssim())
ConstantRootParams(;kwargs...) = Params(assimilation=ConstantNAssim())
ConstantPlant(;kwargs...) = Organism(; params=(ConstantShootParams(), ConstantRootParams()), kwargs...)

    
LeafParams(;kwargs...) = Params(;name=:leaf, translocation=Translocation(destnames=:stem), 
                                assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
StemParams(;kwargs...) = Params(;name=:stem, translocation=Translocation(destnames=(:leaf, :root), proportions=0.5), 
                                assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
RootParams3(;kwargs...) = Params(;name=:root, translocation=Translocation(destnames=:stem), 
                                 assimilation=NAssim(), kwargs...)
LeafVars(;kwargs...) = Vars(;assimilation=CarbonVars(), kwargs...)
StemVars(;kwargs...) = Vars(;assimilation=nothing, kwargs...)
Plant3(;kwargs...) = Organism(;params=(LeafParams(), StemParams(), RootParams()), 
                              vars=(LeafVars(), StemVars(), RootVars()), kwargs...)
