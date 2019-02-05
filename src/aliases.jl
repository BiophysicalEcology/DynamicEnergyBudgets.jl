FLUX = 1.0mol/hr

ShootVars(;kwargs...) = Vars(;assimilation_vars=CarbonVars(), kwargs...)
RootVars(;kwargs...) = Vars(;assimilation_vars=NitrogenVars(), kwargs...)

ShootParamsCNE(;kwargs...) = Params(;name=:shoot, turnover_pars=TurnoverCNE, assimilation_pars=KooijmanSLAPhotosynthesis(), kwargs...)
RootParamsCNE(;kwargs...) = Params(;name=:root, turnover_pars=TurnoverCNE, assimilation_pars=NAssim(), kwargs...)
PlantCNE(;kwargs...) = Plant(;params=(ShootParams(), RootParams()), vars=(ShootVars(), RootVars()), kwargs...)

ShootParamsCN(;kwargs...) = Params(;name=:shoot, assimilation_pars=KooijmanSLAPhotosynthesis(), kwargs...)
RootParamsCN(;kwargs...) = Params(;name=:root, assimilation_pars=NAssim(), kwargs...)
PlantCN(;kwargs...) = Plant(;params=(ShootParamsCN(), RootParamsCN()), vars=(ShootVars(), RootVars()), kwargs...)

FvCBShootParams(;kwargs...) = Params(;name=:shoot, assimilation_pars=FvCBPhotosynthesis(), kwargs...)
FvCBRootParams(;kwargs...) = Params(;name=:root, assimilation_pars=NAssim(), kwargs...)
FvCBLeafParams3(;kwargs...) = Params(;name=:leaf, assimilation_pars=FvCBPhotosynthesis(),
                                    trans_pars=LosslessMultipleTranslocation(destnames=:stem), kwargs...)
FvCBRootParams3(;kwargs...) = Params(;name=:root, trans_pars=LosslessMultipleTranslocation(destnames=:stem), 
                                    assimilation_pars=NAssim(), kwargs...)
FvCBStemParams3(;kwargs...) = Params(;name=:stem, assimilation_pars=nothing, 
                                    trans_pars=LosslessMultipleTranslocation((:leaf,:root), 0.5), kwargs...)
FvCBShootVars(;kwargs...) = Vars(;assimilation_vars=Photosynthesis.PhotoVars(), kwargs...)
FvCBStemVars(;kwargs...) = Vars(;assimilation_vars=nothing, kwargs...)
FvCBLeafVars(;kwargs...) = Vars(;assimilation_vars=Photosynthesis.PhotoVars(), kwargs...)
FvCBRootVars(;kwargs...) = Vars(;assimilation_vars=NitrogenVars(), kwargs...)
FvCBNoRoots(;kwargs...) = Plant(;params=(FvCBShoot()), kwargs...)
FvCBPlant(;kwargs...) = Plant(;params=(FvCBShootParams(), FvCBRootParams()),
                                  vars=(FvCBShootVars(), FvCBRootVars()), kwargs...)
FvCBPlant3(;kwargs...) = Plant(;params=(FvCBLeafParams(), FvCBStemParams(), FvCBRootParams()), 
                              vars=(FvCBLeafVars(), FvCBStemVars(), FvCBRootVars()), kwargs...)

ConstantShootParams(;kwargs...) = Params(turnover_pars=TurnoverCNE, assimilation_pars=ConstantCAssim())
ConstantRootParams(;kwargs...) = Params(turnover_pars=TurnoverCNE, assimilation_pars=ConstantNAssim())
ConstantPlantCNE(;kwargs...) = Plant(; params=(ConstantShootParamsCNE(), ConstantRootParamsCNE()), kwargs...)

ConstantShootParamsCN(;kwargs...) = Params(assimilation_pars=ConstantCAssim())
ConstantRootParamsCN(;kwargs...) = Params(assimilation_pars=ConstantNAssim())
ConstantPlantCN(;kwargs...) = Plant(; params=(ConstantShootParamsCN(), ConstantRootParamsCN()), kwargs...)

    
LeafParams(;kwargs...) = Params(;name=:leaf, trans_pars=LosslessMultipleTranslocation(destnames=:stem), 
                                assimilation_pars=KooijmanSLAPhotosynthesis(), kwargs...)
StemParams(;kwargs...) = Params(;name=:stem, trans_pars=LosslessMultipleTranslocation(destnames=(:leaf, :root), proportions=0.5), 
                                assimilation_pars=nothing, kwargs...)
RootParams3(;kwargs...) = Params(;name=:root, trans_pars=LosslessMultipleTranslocation(destnames=:stem), 
                                 assimilation_pars=NAssim(), kwargs...)
LeafVars(;kwargs...) = Vars(;assimilation_vars=CarbonVars(), kwargs...)
StemVars(;kwargs...) = Vars(;assimilation_vars=nothing, kwargs...)
Plant3(;kwargs...) = Plant(;params=(LeafParams(), StemParams(), RootParams()), 
                              vars=(LeafVars(), StemVars(), RootVars()), kwargs...)
