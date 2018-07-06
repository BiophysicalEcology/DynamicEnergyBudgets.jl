Shoot(;kwargs...) = Organ(params=Params(assimilation=KooijmanSLAPhotosynthesis()), 
                              vars=Vars(assimilation=nothing), kwargs...)
Root(;kwargs...) = Organ(params=Params(assimilation=N_Assimilation()), 
                             vars=Vars(assimilation=nothing), kwargs...)
Plant(;kwargs...) = Organism(organs=(Shoot(), Root))


ConstantShoot(;kwargs...) = Organ(params=Params(assimilation=ConstantCarbonAssimilation()), 
                        vars=Vars(assimilation=nothing), kwargs...)
ConstantRoot(;kwargs...) = Organ(params=Params(assimilation=ConstantNitrogenAssimilation()), 
                       vars=Vars(assimilation=nothing), kwargs...)
ConstantPlant(;kwargs...) = Organism(; organs=(ConstantShoot(), ConstantRoot()), kwargs...) 


UntypedVars(;kwargs...) = default_kw(Vars{Any,Any,Any,Any,Any}; kwargs...)
UntypedOrgan(state::S, name::N, params::P, shared::Sh, vars::V) where {S,N,P,Sh,V} = begin
    # Get flux units from params, instead of explicitly.
    one_flux = oneunit_flux(params, state)
    J = build_J(one_flux, state; typ = Any)
    J1 = build_J1(one_flux, state; typ = Any)
    Organ{S,N,P,Sh,V,typeof(J),typeof(J1)}(state, name, params, shared, vars, J, J1)
end
UntypedShoot(;kwargs...) = Organ(params=Params(assimilation=ConstantCarbonAssimilation()), 
                                 vars=UntypedVars(assimilation=nothing), kwargs...)
UntypedRoot(;kwargs...) = Organ(params=Params(assimilation=ConstantNitrogenAssimilation()), 
                                vars=UntypedVars(assimilation=nothing), kwargs...)
UntypedPlant(;kwargs...) = Organism(time=0:1:1000, organs=(UntypedShoot(), UntypedRoot()), kwargs...) 


FvCBShootParams(;kwargs...) = Params(; assimilation=PhotoParams(), kwargs...)
FvCBRootParams(;kwargs...) = Params(; assimilation=N_Assimilation(), kwargs...)
FvCBShootVars(;kwargs...) = Vars(; assimilation=PhotoVars(), kwargs...)
FvCBRootVars(;kwargs...) = Vars(; assimilation=NitrogenVars(), kwargs...)
FvCBShoot(;kwargs...) = Organ(; params=FvCBShootParams, vars=FvCBRootVars(), kwargs...)
FvCBRoot(;kwargs...) = Organ(; params=FvCBRootParams, vars=FvCBRootVars(), kwargs...)
FvCBPlant(;kwargs...) = Organism(; organs=(FvCBShoot(), FvCBRoot()), kwargs...) 
