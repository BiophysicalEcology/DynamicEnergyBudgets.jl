FLUX = 1.0mol/hr

ShootParams(;kwargs...) = ParamsCNE(;name=:shoot, assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
RootParams(;kwargs...) = ParamsCNE(;name=:root, assimilation=NAssim(), kwargs...)
ShootVars(;kwargs...) = Vars(;assimilation=CarbonVars(), kwargs...)
RootVars(;kwargs...) = Vars(;assimilation=NitrogenVars(), kwargs...)
Records(params; vars=(ShootVars(), RootVars()), t=1:1:1, val=FLUX) = begin (Records(params[1], vars[1], t, val, typeof(val)), 
     Records(params[2], vars[2], t, val, typeof(val)))
end
Plant(;kwargs...) = Organism(;params=(ShootParams(), RootParams()), vars=(ShootVars(), RootVars()), kwargs...)

ShootParamsCN(;kwargs...) = ParamsCN(;name=:shoot, assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
RootParamsCN(;kwargs...) = ParamsCN(;name=:root, assimilation=NAssim(), kwargs...)
PlantCN(;kwargs...) = Organism(;params=(ShootParamsCN(), RootParamsCN()), vars=(ShootVars(), RootVars()), kwargs...)

NoMaturityShootParams(;kwargs...) = ShootParams(;maturity=nothing, kwargs...)
NoMaturityRootParams(;kwargs...) = RootParams(;maturity=nothing, kwargs...)
NoMaturityPlant(;kwargs...) = Organism(;params=(NoMaturityShootParams(), NoMaturityRootParams()), vars=(ShootVars(), RootVars()), kwargs...)


FvCBShootParams(;kwargs...) = ParamsCNE(;name=:shoot, assimilation=FvCBPhotosynthesis(), kwargs...)
FvCBRootParams(;kwargs...) = ParamsCNE(;name=:root, assimilation=NAssim(), kwargs...)
FvCBLeafParams3(;kwargs...) = ParamsCNE(;name=:leaf, assimilation=FvCBPhotosynthesis(),
                                    translocation=LosslessMultipleTranslocation(destnames=:stem), kwargs...)
FvCBRootParams3(;kwargs...) = ParamsCNE(;name=:root, translocation=LosslessMultipleTranslocation(destnames=:stem), 
                                    assimilation=NAssim(), kwargs...)
FvCBStemParams3(;kwargs...) = ParamsCNE(;name=:stem, assimilation=nothing, 
                                    translocation=LosslessMultipleTranslocation((:leaf,:root), 0.5), kwargs...)
FvCBShootVars(;kwargs...) = Vars(;assimilation=Photosynthesis.PhotoVars(), kwargs...)
FvCBStemVars(;kwargs...) = Vars(;assimilation=nothing, kwargs...)
FvCBLeafVars(;kwargs...) = Vars(;assimilation=Photosynthesis.PhotoVars(), kwargs...)
FvCBRootVars(;kwargs...) = Vars(;assimilation=NitrogenVars(), kwargs...)
FvCBNoRoots(;kwargs...) = Organism(;params=(FvCBShoot()), kwargs...)
FvCBPlant(;kwargs...) = Organism(;params=(FvCBShootParams(), FvCBRootParams()),
                                  vars=(FvCBShootVars(), FvCBRootVars()), kwargs...)
FvCBPlant3(;kwargs...) = Organism(;params=(FvCBLeafParams(), FvCBStemParams(), FvCBRootParams()), 
                              vars=(FvCBLeafVars(), FvCBStemVars(), FvCBRootVars()), kwargs...)


ConstantShootParamsCNE(;kwargs...) = ParamsCNE(assimilation=ConstantCAssim())
ConstantRootParamsCNE(;kwargs...) = ParamsCNE(assimilation=ConstantNAssim())
ConstantPlantCNE(;kwargs...) = Organism(; params=(ConstantShootParamsCNE(), ConstantRootParamsCNE()), kwargs...)

ConstantShootParamsCN(;kwargs...) = ParamsCN(assimilation=ConstantCAssim())
ConstantRootParamsCN(;kwargs...) = ParamsCN(assimilation=ConstantNAssim())
ConstantPlantCN(;kwargs...) = Organism(; params=(ConstantShootParamsCN(), ConstantRootParamsCN()), kwargs...)

    
LeafParams(;kwargs...) = ParamsCNE(;name=:leaf, translocation=LosslessMultipleTranslocation(destnames=:stem), 
                                assimilation=KooijmanSLAPhotosynthesis(), kwargs...)
StemParams(;kwargs...) = ParamsCNE(;name=:stem, translocation=LosslessMultipleTranslocation(destnames=(:leaf, :root), proportions=0.5), 
                                assimilation=nothing, kwargs...)
RootParams3(;kwargs...) = ParamsCNE(;name=:root, translocation=LosslessMultipleTranslocation(destnames=:stem), 
                                 assimilation=NAssim(), kwargs...)
LeafVars(;kwargs...) = Vars(;assimilation=CarbonVars(), kwargs...)
StemVars(;kwargs...) = Vars(;assimilation=nothing, kwargs...)
Plant3(;kwargs...) = Organism(;params=(LeafParams(), StemParams(), RootParams()), 
                              vars=(LeafVars(), StemVars(), RootVars()), kwargs...)
