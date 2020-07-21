
"""
Astract supertype for organs parameters. 
Extend to change the components that are specific to each organ.
"""
abstract type AbstractParams end

"""
Model parameters that vary between organs
"""
@default_kw @flattenable @selectable struct Params{As,Sh,Al,Ma,Tr,Re,Ge,Pr} <: AbstractParams
    # Field               | Default                    | _     | Selectable Types
    assimilation_pars::As | ConstantCAssim()           | _     | Union{Nothing,AbstractAssim}
    scaling_pars::Sh      | nothing                    | _     | Union{Nothing,AbstractScaling}
    allometry_pars::Al    | nothing                    | _     | Union{Nothing,AbstractAllometry}
    maturity_pars::Ma     | nothing                    | _     | Union{Nothing,AbstractMaturity}
    activetrans_pars::Tr  | nothing                    | _     | Union{Nothing,ActiveTranslocation}
    passivetrans_pars::Re | LosslessPassiveTranslocation() | _ | PassiveTranslocation
    germination_pars::Ge  | nothing                    | _     | Union{Nothing,AbstractGermination}
    production_pars::Pr   | nothing                    | _     | Union{Nothing,AbstractProduction}
end

for fn in fieldnames(Params)
    @eval $fn(p::Params) = p.$fn
end


"""
Astract supertype for shared parameters. 
Extend to change the components that are shared.
"""
abstract type AbstractSharedParams end

"""
    SharedParams(su_pars, core_pars, resorption_pars, tempcorr_pars, catabolism_pars)

Model parameters shared between organs.

# FieldMetadata macros
- `@default_kw` provides default values and constructors, 
- `@selectable` provides the set of types that can be used for the parameter,
  so that they can be selected in a live interface.
"""
@default_kw @selectable struct SharedParams{SU,Co,Fe,Te,Ca} <: AbstractSharedParams
    # Field               | Default                   | Selectable Types
    su_pars::SU           | ParallelComplementarySU() | AbstractSynthesizingUnit
    core_pars::Co         | DEBCore()                 | _
    resorption_pars::Fe   | nothing                   | Union{Nothing,AbstractResorption}
    tempcorr_pars::Te     | nothing                   | Union{Nothing,AbstractTemperatureCorrection}
    catabolism_pars::Ca   | CatabolismCN()            | AbstractCatabolism
end

for fn in fieldnames(SharedParams)
    @eval $fn(p::SharedParams) = p.$fn
end


###########################################################################################
# Variables

"""
Model variables. 
Allow storing and accessing variables for use by multiple components.
"""
abstract type AbstractVars end

"""
    PlottableVars()

Plottable model variables. These are vectors witih values for each time-step,
to allow plotting and model introspection.
"""
@udefault_kw @units @plottable struct PlottableVars{F,R,E,T,WP,H,TS} <: AbstractVars
    scaling::F        | [0.0] | _                 | _
    rate::R           | [0.0] | mol*mol^-1*d^-1   | _
    E_ctb::E          | [0.0] | mol/hr            | _
    θE::F             | [0.0] | _                 | _
    temp::T           | [0.0] | K                 | _
    tempcorrection::F | [0.0] | _                 | _
    swp::WP           | [0.0] | kPa               | _
    soilcorrection::F | [0.0] | _                 | _
    height::H         | [0.0] | m                 | _
    tstep::TS         | [1]   | _                 | false
end


"""
    Vars()

Mutable struct to allow storing variables
for use by multiple components.
"""
@udefault_kw @units mutable struct Vars{F,R,E,T,WP,H} <: AbstractVars
    scaling::F        | 0.0   | _
    rate::R           | 0.0   | mol*mol^-1*d^-1
    θE::F             | 0.0   | _
    E_ctb::E          | 0.0   | mol/hr          
    temp::T           | 0.0   | K
    tempcorrection::F | 0.0   | _
    swp::WP           | 0.0   | kPa
    soilcorrection::F | 0.0   | _
    height::H         | 0.0   | m
end

tstep(v::PlottableVars) = v.tstep[1]
tstep(v::Vars) = nothing

set_tstep!(v::PlottableVars, val) = v.tstep[1] = val
set_tstep!(v::Vars, val) = nothing

depth(v) = height(v)

build_vars(vars::Vars, tspan) = vars
build_vars(vars::PlottableVars, tspan) = begin
    len = length(tspan)
    len <= length(vars.rate) && return vars

    for fname in fieldnames(typeof(vars))
        varvec = getfield(vars, fname)
        append!(varvec, fill(getfield(vars, fname)[1], len - length(varvec)))
    end
    vars
end

###########################################################################################
# Organs and Organisms

"""
Abstract supertype for organs. Inherit from it if you need to difine
behaviour diferent to that or [`Organ`](@ref).
"""
abstract type AbstractOrgan{P,S} end

params(o::AbstractOrgan) = o.params
shared(o::AbstractOrgan) = o.shared
vars(o::AbstractOrgan) = o.vars
flux(o::AbstractOrgan) = o.J

#= κ functions allow modular addition of destinations for catabolised flux, 
by default maturity and active translocation. If those components are
not used the flux fraction will be allocated to κsoma. =#

κtra(o::AbstractOrgan) = κtra(activetrans_pars(o))
κtra(o::Nothing) = 0.0

κmat(o::AbstractOrgan) = κmat(maturity_pars(o))
κmat(::Nothing) = 0.0

κsoma(o::AbstractOrgan) = oneunit(κtra(o)) - κtra(o) - κmat(o)

# Define `scaling` and `setscaling` etc. methods
for fn in fieldnames(Vars)
    setfn = Symbol.(:set_, fn, :!)
    @eval $fn(o::AbstractOrgan) = $fn(vars(o))
    @eval $setfn(o::AbstractOrgan, x) = $setfn(vars(o), x)
    @eval @inline ($fn)(vars::PlottableVars) = vars.$fn[tstep(vars)]
    @eval @inline ($setfn)(vars::PlottableVars, val) = vars.$fn[tstep(vars)] = val
    @eval @inline ($fn)(vars::Vars) = vars.$fn
    @eval @inline ($setfn)(vars::Vars, val) = vars.$fn = val
end

#= Forward variable and parameter getter/setter methods
so they can be acessesed/set directly from the organ,
without knowning where they are actually stored. =#
tstep(o::AbstractOrgan) = tstep(vars(o))
set_tstep!(o::AbstractOrgan, t) = set_tstep!(vars(o), t)

for fn in fieldnames(DEBCore)
    @eval $fn(p::SharedParams) = $fn(core_pars(p))
    @eval $fn(o::AbstractOrgan{<:Any,<:SharedParams}) = $fn(shared(o))
end

for fn in fieldnames(SharedParams)
    @eval $fn(o::AbstractOrgan{<:Any,<:SharedParams}) = $fn(shared(o))
end

for fn in fieldnames(Params)
    @eval $fn(o::AbstractOrgan{<:Params}) = $fn(params(o))
end

"""
    Organ(params, shared, vars, J)

Basic model components. For a plants, organs might be roots, stem and leaves
"""
struct Organ{P,S,V,F} <: AbstractOrgan{P,S}
    params::P
    shared::S
    vars::V
    J::F
end
"""
    Organ(params, shared, records)

Construct an organ from parameters, shared parameters and
views into records arrays for vaiable and flux matrices.
"""
Organ(params::AbstractParams, shared::AbstractSharedParams, records) = begin
    vars = records.vars
    set_tstep!(vars, 1)
    J = view(records.J, Ti(1))
    Organ(params, shared, vars, J)
end


abstract type AbstractRecords end

"""
    Records(vars, J)

Time series of mutable variables and flux for ploting and analysis

These are sliced with `view` for each timestep. An effecient implementation
may use a single view repeatedly, losing ability to plot values over time.
"""
@plottable struct PlottableRecords{V,F} <: AbstractRecords
    vars::V | true
    J::F    | false
end
PlottableRecords(vars::PlottableVars, J::AbstractArray) =
    PlottableRecords{map(typeof,(vars,J))...}(vars, J)
PlottableRecords(vars::PlottableVars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) = begin
    vars = build_vars(vars, tspan)
    J = build_flux(fluxval, fluxaxes..., tspan)
    PlottableRecords(vars, J)
end

@plottable struct Records{V,F} <: AbstractRecords
    vars::V | false
    J::F    | false
end
Records(vars::Vars, J::AbstractArray) =
    Records{map(typeof,(vars,J))...}(vars, J)
Records(vars::Vars, fluxval::Number, fluxaxes::Tuple) = begin
    J = build_flux(fluxval, fluxaxes...)
    Records{map(typeof, (vars,J))...}(vars, J)
end

build_flux(fluxval, x::Tuple, y::Tuple) = begin
    dims = X(Val(x)), Y(Val(y))
    A = zeros(typeof(fluxval), map(length, dims)...)
    DimensionalArray(A, dims)
end
build_flux(fluxval, x::Tuple, y::Tuple, time::AbstractRange) = begin
    dims = X(Val(x)), Y(Val(y)), Ti(time)
    A = zeros(typeof(fluxval), map(length, dims)...)
    DimensionalArray(A, dims)
end



abstract type AbstractOrganism end

params(o::AbstractOrganism) = o.params
shared(o::AbstractOrganism) = o.shared
records(o::AbstractOrganism) = o.records
organs(o::AbstractOrganism) = o.organs
environment(o::AbstractOrganism) = o.environment
environment_start(o::AbstractOrganism) = o.environment_start
dead(o::AbstractOrganism) = o.dead[]
set_dead!(o::AbstractOrganism, val) = o.dead[] = val

"""
    define_organs(o::AbstractOrganism, t)

Organs are constructed with views of Records and J Arrays at time t
"""
define_organs(o::AbstractOrganism, t) =
    define_organs(params(o), shared(o), records(o), t)
define_organs(params::Tuple, shared, records::Tuple, t) =
    map((p, r) -> Organ(p, shared, r), params, records)

update_organs(organs, t) = update_organs(vars(first(organs)), organs, t)
update_organs(::Vars, organs, t) = organs
update_organs(::PlottableVars, organs, t) = begin
    i = floor(Int, ustrip(t))
    map(organs) do organ
        organ.vars.tstep[1] = i
        J = organ.J
        data = J.data
        @set! data.indices[3] = i
        @set! data.offset1 = (i - 1) * length(data.indices[1]) * length(data.indices[2]) 
        # Setfield.jl isn't type-stable, DD `rebuild` is better
        @set! organ.J = DimensionalData.rebuild(J, data)
        organ
    end
end

"""
    Plant(params, shared, records, environment, environment_start, dead)

Basic plant model parameters.
"""
@flattenable @description mutable struct Plant{P,S,R,O,E,ES,D} <: AbstractOrganism
    params::P             | true  | "Model parameters"
    shared::S             | true  | "Parameters shared between organs"
    records::R            | false | "Plotable variables stored in arrays and sliced for each timestep on demand"
    organs::O             | false | ""
    environment::E        | false | "Environment object, provides environmental variables for each timestep"
    environment_start::ES | false | "Start index of environmental data"
    dead::D               | false | "`Bool` flag: has the plant died"
    Plant(params::P, shared::S, records::R, organs::O, environment::E, environment_start::ES, dead::D) where {P,S,R,O,E,ES,D} = begin
        organs = define_organs(params, shared, records, 0)
        new{P,S,R,typeof(organs),E,ES,D}(params, shared, records, organs, environment, environment_start, dead)
    end
end

"""
    Plant(states=(:V, :C, :N),
          transformations=(:asi, :gro, :mai, :rej, :res),
          catstates=(:CN, :C, :N, :E),
          cattransformations=(:ctb,),
          params=(ShootParamsCN(), RootParamsCN()),
          vars=(Vars(), Vars()),
          shared=SharedParams(),
          records=nothing,
          environment=nothing,
          time=0.0hr:1.0hr:8760.0hr,
          environment_start=Ref(1.0hr),
          dead=Ref(false))

Outer construtor for defaults
"""
Plant(; states=(:V, :C, :N),
        transformations=(:asi, :gro, :mai, :rej, :tra, :res),
        params=(
            Params(assimilation_pars=KooijmanSLAPhotosynthesis()),
            Params(assimilation_pars=NAssim()),
        ),
        vars=(Vars(), Vars()),
        shared=SharedParams(),
        records=nothing,
        environment=nothing,
        time=0.0hr:1.0hr:8760.0hr,
        environment_start=Ref(0.0hr),
        dead=Ref(false)
      ) = begin
    fluxaxes = states, transformations
    fluxval = 1.0mol/hr
    records = build_records(records, vars, fluxval, fluxaxes, time)
    Plant(params, shared, records, nothing, environment, environment_start, dead)
end

build_records(records::AbstractRecords, args...) = records
build_records(records::Nothing, vars::Tuple, fluxval, fluxaxes, tspan) =
    map(vars) do v
        build_records(v, fluxval, fluxaxes, tspan)
    end
build_records(vars::Vars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) =
    Records(vars, fluxval, fluxaxes)
build_records(vars::PlottableVars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) =
    PlottableRecords(vars, fluxval, fluxaxes, tspan)



# Define a SubArray constructor so we can increment the time index
ConstructionBase.constructorof(::Type{A}) where A<:SubArray{T,N,P,I,L} where {T,N,P,I,L} =
    (parent, indices, offset1, stride1) -> begin
        SubArray{eltype(parent),ndims(parent),typeof(parent),typeof(indices),L}(parent, indices, offset1, stride1)
    end
