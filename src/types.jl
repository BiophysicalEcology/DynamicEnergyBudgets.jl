

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

" Model parameters that vary between organs "
@selectable @flattenable @default_kw struct Params{As,Sh,Al,Ma,Tr,Re,Ge,Pr} <: AbstractParams
    # Field               | Default                | _     | Selectable Types
    name::Symbol          | :organ                 | false | _
    rate_formula          | FZeroRate()            | _     | _
    assimilation_pars::As | ConstantCAssim()       | _     | Union{Nothing,AbstractAssim}
    shape_pars::Sh        | Plantmorph()           | _     | Union{Nothing,AbstractShape}
    allometry_pars::Al    | nothing                | _     | Union{Nothing,AbstractAllometry}
    maturity_pars::Ma     | nothing                | _     | Union{Nothing,AbstractMaturity}
    trans_pars::Tr        | nothing                | _     | Union{Nothing,AbstractTranslocation}
    rejection_pars::Re    | LosslessRejection()    | _     | AbstractRejection
    germination_pars::Ge  | ThresholdGermination() | _     | Union{Nothing,AbstractGermination}
    production_pars::Pr   | Production()           | _     | Union{Nothing,AbstractProduction}
end


rate_formula(p) = p.rate_formula
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


" Model parameters shared between organs "
@selectable @udefault_kw struct SharedParams{SU,Co,Fe,Te,Ca,Mt}
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


" Model variables "
@plottable @units @udefault_kw struct Vars{F,MoMoD,C,WP,M,T}
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
    set = Symbol.(:set_, f, :!)
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

" Basic model components. For a plants, organs might be roots, stem and leaves "
struct Organ{P,S,V,F,F1} <: AbstractOrgan
    params::P
    shared::S
    vars::V
    J::F
    J1::F1
end

" update vars and flux to record for time t "
define_organs(o, t) = define_organs(params(o), shared(o), records(o), t)
define_organs(params::Tuple{P,Vararg}, shared, records::Tuple{R,Vararg}, t) where {P,R} = begin
    rec = records[1]
    t = calc_tstep(rec.vars, t)
    rec.vars.t[1] = t
    vJ = view(rec.J, :, :, t)
    vJ1 = view(rec.J1, :, :, t)
    J = LArray{Tuple{STATE,TRANS}}(vJ)
    J1 = LArray{Tuple{STATE1,TRANS1}}(vJ1)

    organ = Organ(params[1], shared, rec.vars, J, J1)
    (organ, define_organs(tail(params), shared, tail(records), t)...)
end
define_organs(params::Tuple{}, shared, records::Tuple{}, t) = ()


"Records of mutable variables and flux for ploting and analysis"
@plottable struct Records{V,F,F1}
    vars::V | true
    J::F    | true
    J1::F1  | true 
end
"Constructor for records. Arrays use the length of the current timespan"
Records(params, vars, time, val, typ) = begin
    vars = build_vars(vars, time)
    J = build_J(val, typ, time)
    J1 = build_J1(val, typ, time)
    Records(vars, J, J1)
end

build_J(one_flux, T, time) = build_flux(one_flux, T, STATE, TRANS, time)

build_J1(one_flux, T, time) = build_flux(one_flux, T, STATE1, TRANS1, time)

build_flux(one_flux, T, x, y, time) = zeros(typeof(one_flux), length(x), length(y), length(time))


abstract type AbstractOrganism end

dead(o::AbstractOrganism) = o.dead[]
set_dead!(o::AbstractOrganism, val) = o.dead[] = val
environment(o::AbstractOrganism) = o.environment

vars(o::AbstractOrgan) = o.vars
flux(o::AbstractOrgan) = o.J
flux1(o::AbstractOrgan) = o.J1

"An organism, made up of organs"
@flattenable mutable struct Plant{P,S,R,E,ES,D} <: AbstractOrganism
    params::P             | true
    shared::S             | true
    records::R            | false
    environment::E        | true
    environment_start::ES | false
    dead::D               | false
end

"Outer construtor for defaults"
Plant(; params = (ShootParamsCN(), RootParamsCN()),
        vars = (Vars(), Vars()),
        shared = SharedParams(),
        records = nothing,
        environment = nothing,
        time = 0hr:1hr:8760hr,
        environment_start = Ref(1hr),
        dead = Ref(false)
      ) = begin
    if records == nothing
        records = []
        for i = 1:length(params)
            push!(records, Records(params[i], vars[i], time, 1.0mol/hr, typeof(1.0mol/hr)))
        end
        records = (records...,)
    end
    Plant(params, shared, records, environment, environment_start, dead)
end
