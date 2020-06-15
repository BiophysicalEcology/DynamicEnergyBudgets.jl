

# So we don't have to depend on all of Lazy.jl
macro forward(ex, fs)
  @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
  T = esc(T)
  fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
  :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
       for f in fs]...);
    nothing)
end


abstract type AbstractParams end

"""
Model parameters that vary between organs
"""
@default_kw @flattenable @selectable struct Params{As,Sh,Al,Ma,Tr,Re,Ge,Pr} <: AbstractParams
    # Field               | Default             | _     | Selectable Types
    name::Symbol          | :organ              | false | _
    assimilation_pars::As | ConstantCAssim()    | _     | Union{Nothing,AbstractAssim}
    shape_pars::Sh        | Plantmorph()        | _     | Union{Nothing,AbstractShape}
    allometry_pars::Al    | nothing             | _     | Union{Nothing,AbstractAllometry}
    maturity_pars::Ma     | nothing             | _     | Union{Nothing,AbstractMaturity}
    trans_pars::Tr        | nothing             | _     | Union{Nothing,AbstractTranslocation}
    rejection_pars::Re    | LosslessRejection() | _     | AbstractRejection
    germination_pars::Ge  | nothing             | _     | Union{Nothing,AbstractGermination}
    production_pars::Pr   | nothing             | _     | Union{Nothing,AbstractProduction}
end

assimilation_pars(p) = p.assimilation_pars
shape_pars(p) = p.shape_pars
allometry_pars(p) = p.allometry_pars
maturity_pars(p) = p.maturity_pars
trans_pars(p) = p.trans_pars
rejection_pars(p) = p.rejection_pars
germination_pars(p) = p.germination_pars
production_pars(p) = p.production_pars

# turnover_pars(p) = p.turnover_pars

n_N_P(p) = production_pars(p).n_N_P
y_V_E(p) = core_pars(p).y_V_E
y_E_EC(p) = core_pars(p).y_E_EC
y_E_EN(p) = core_pars(p).y_E_EN
n_N_V(p) = core_pars(p).n_N_V
n_N_E(p) = core_pars(p).n_N_E
n_N_EC(p) = core_pars(p).n_N_EC
n_N_EN(p) = core_pars(p).n_N_EN
w_V(p) = core_pars(p).w_V
w_C(p) = core_pars(p).w_C
w_N(p) = core_pars(p).w_N
w_E(p) = core_pars(p).w_E

#    %   W   Rel n atoms
# C  45  12  30857
# H  6   1   58600
# O  45  16  27000
# N  2   14   1028
# K  1   39    246
# X  1   40    300ish

abstract type AbstractSharedParams end

"""
    SharedParams(su_pars, core_pars, resorption_pars, tempcorr_pars, catabolism_pars, maintenance_pars)

Model parameters shared between organs.
"""
@udefault_kw @selectable struct SharedParams{SU,Co,Fe,Te,Ca,Mt} <: AbstractSharedParams
    # Field               | Default                   | Selectable Types
    su_pars::SU           | ParallelComplementarySU() | AbstractSynthesizingUnit
    core_pars::Co         | DEBCore()                 | _
    resorption_pars::Fe   | nothing                   | Union{Nothing,AbstractResorption}
    tempcorr_pars::Te     | nothing                   | Union{Nothing,AbstractTemperatureCorrection}
    catabolism_pars::Ca   | CatabolismCN()            | AbstractCatabolism
    maintenance_pars::Mt  | Maintenance()             | AbstractMaintenance
end

su_pars(p) = p.su_pars
core_pars(p) = p.core_pars
resorption_pars(p) = p.resorption_pars
tempcorr_pars(p) = p.tempcorr_pars
catabolism_pars(p) = p.catabolism_pars
maintenance_pars(p) = p.maintenance_pars

###########################################################################################
# Variables

abstract type AbstractVars end

"""
Model variables
"""
@udefault_kw @units @plottable struct Vars{F,MoMoD,C,WP,M,T} <: AbstractVars
    shape::F             | [0.0]   | _                 | _
    rate::MoMoD          | [0.0]   | mol*mol^-1*d^-1   | _
    θE::F                | [0.0]   | _                 | _
    temp::C              | [25.0]  | K                 | _
    tempcorrection::F    | [1.0]   | _                 | _
    swp::WP              | [1.0]   | kPa               | _
    soilcorrection::F    | [1.0]   | _                 | _
    height::M            | [0.0]   | m                 | _
    t::T                 | [1]     | _                 | false
end

# Define `shape` and `setshape` etc. methods
for field in [:shape, :rate, :temp, :θE, :tempcorrection, :soilcorrection, :height, :swp]
    set = Symbol.(:set_, field, :!)
    @eval @inline ($field)(vars) = vars.$field[tstep(vars)]
    @eval @inline ($set)(vars, val) = vars.$field[tstep(vars)] = val
end


tstep(v) = v.t[1]
set_tstep!(v, val) = v.t[1] = val
depth(v) = height(v)

build_vars(vars, time) = begin
    len = length(time)
    len == length(vars.rate) && return vars

    fields = []
    for fname in fieldnames(typeof(vars))
        ft = fieldtype(typeof(vars), fname)
        if ft <: AbstractArray
            push!(fields, fill(getfield(vars, fname)[1], len))
        else
            push!(fields, getfield(vars, fname))
        end
    end
    typeof(vars).name.wrapper(fields...)
end

###########################################################################################
# Organs and Organisms

abstract type AbstractOrgan end

params(o::AbstractOrgan) = o.params
shared(o::AbstractOrgan) = o.shared
vars(o::AbstractOrgan) = o.vars
flux(o::AbstractOrgan) = o.J
flux1(o::AbstractOrgan) = o.J1

j_E_mai(o::AbstractOrgan) = j_E_mai(maintenance_pars(o))

κtra(o::AbstractOrgan) = κtra(trans_pars(o))
κtra(o::Nothing) = 0.0

κmat(o::AbstractOrgan) = κmat(maturity_pars(o))
κmat(::Nothing) = 0.0

κsoma(o::AbstractOrgan) = oneunit(κtra(o)) - κtra(o) - κmat(o)

@forward AbstractOrgan.vars θE, temp, set_temp!, set_swp!, swp, set_soilcorrection!, soilcorrection, tempcorrection, set_tempcorrection!,
         height, set_height!, rate, set_rate!, shape, set_shape!, tstep, set_tstep!

@forward AbstractOrgan.params rate_formula, assimilation_pars, shape_pars, allometry_pars, maturity_pars,
                              trans_pars, production_pars, rejection_pars, germination_pars, turnover_pars

@forward AbstractOrgan.shared maintenance_pars, resorption_pars, su_pars, tempcorr_pars, catabolism_pars, core_pars,
                              y_V_E, y_E_EC, y_E_EN, n_N_P, n_N_V, n_N_E, n_N_EC, n_N_EN, w_V, w_C, w_N, w_E

"""
    Organ(params, shared, vars, J, J1)

Basic model components. For a plants, organs might be roots, stem and leaves
"""
struct Organ{P,S,V,F,F1} <: AbstractOrgan
    params::P
    shared::S
    vars::V
    J::F
    J1::F1
end
"""
    Organ(params, shared, records, t)

Construct an organ from parameters, shared parameters and
views into records arrays for vaiable and flux matrices.
"""
Organ(params::AbstractParams, shared::AbstractSharedParams, records, t) = begin
    vars = records.vars
    t = calc_tstep(vars, t)
    vars.t[1] = t
    J = view(records.J, Ti(t))
    J1 = view(records.J1, Ti(t))
    Organ(params, shared, vars, J, J1)
end


"""
    Records(vars, J, J1)

Time series of mutable variables and flux for ploting and analysis

These are sliced with `view` for each timestep. An effecient implementation
may use a single view repeatedly, losing ability to plot values over time.
"""
@plottable struct Records{V,F,F1}
    vars::V | true
    J::F    | true
    J1::F1  | true
end
"""
    Records(params::AbstractParams, vars, time, fluxval)

Constructor for records. Arrays use the length of the current timespan.
"""
Records(states, trans, catstates, cattrans, vars, time, fluxval) = begin
    vars = build_vars(vars, time)
    J = build_flux(fluxval, states, trans, time)
    J1 = build_flux(fluxval, catstates, cattrans, time)
    Records(vars, J, J1)
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
environment(o::AbstractOrganism) = o.environment
environment_start(o::AbstractOrganism) = o.environment_start
dead(o::AbstractOrganism) = o.dead[]
set_dead!(o::AbstractOrganism, val) = o.dead[] = val

"""
    define_organs(o::AbstractOrganism, t)

Organs are constructed with views of Records and J/J1 Arrays at time t
"""
define_organs(o::AbstractOrganism, t) =
    map((p, r) -> Organ(p, shared(o), r, t), params(o), records(o))


"""
    Plant(params, shared, records, environment, environment_start, dead)

Plant model parameters.
"""
@flattenable @description mutable struct Plant{P,S,R,E,ES,D} <: AbstractOrganism
    params::P             | true  | "Model parameters"
    shared::S             | true  | "Parameters shared between organs"
    records::R            | false | "Plotable bariables stored in arrays and sliced for each timestep on demand"
    environment::E        | false | "Environment object, provides environmental variables for each timestep"
    environment_start::ES | false | "Start index of environmental data"
    dead::D               | false | "`Bool` flag: has the plant died"
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
        catstates=(:CN, :C, :N, :E),
        cattransformations=(:ctb,),
        params=(ShootParamsCN(), RootParamsCN()),
        vars=(Vars(), Vars()),
        shared=SharedParams(),
        records=nothing,
        environment=nothing,
        time=0.0hr:1.0hr:8760.0hr,
        environment_start=Ref(1.0hr),
        dead=Ref(false)
      ) = begin
    if records == nothing
        records = []
        for i = 1:length(params)
            push!(records, Records(states, transformations, catstates, cattransformations, vars[i], time, 1.0mol/hr))
        end
        records = (records...,)
    end
    Plant(params, shared, records, environment, environment_start, dead)
end
