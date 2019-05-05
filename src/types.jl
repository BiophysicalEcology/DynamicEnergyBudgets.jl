abstract type AbstractDEBCore end

@columns struct DEBCore{MoMo,GMo} <: AbstractDEBCore
    y_V_E::MoMo  | 0.7   | _        | Beta(2.0, 2.0)  | [0.0, 1.0]   | _     | "From reserve to structure"
    y_E_EC::MoMo | 0.7   | _        | Gamma(2.0, 2.0) | [1e-6, 1.0]  | false | "From C-reserve to reserve, using nitrate"
    y_E_EN::MoMo | 30.0  | _        | Gamma(2.0, 2.0) | [1.0, 50.0]  | false | "From N-reserve to reserve"
    n_N_V::MoMo  | 0.03  | _        | Gamma(2.0, 2.0) | [0.0, 0.1]   | _     | "N/C in structure"
    n_N_E::MoMo  | 0.025 | _        | Gamma(2.0, 2.0) | [0.0, 0.1]   | _     | "N/C in reserve"
    w_V::GMo     | 25.0  | g*mol^-1 | Gamma(2.0, 2.0) | [15.0, 40.0] | _     | "Mol-weight of shoot structure"
    # w_N::GMo     | 25.0  | g*mol^-1 | Gamma(2.0, 2.0) | [15.0, 40.0] | _     | "Mol-weight of shoot N-reserve"
    # w_C::GMo   | 25.0  | g*mol^-1 | Gamma(2.0, 2.0) | [12.0, 40.0] | _     | "Mol-weight of shoot C-reserve"
    # w_E::GMo   | 25.0  | g*mol^-1 | Gamma(2.0, 2.0) | [15.0, 40.0] | _     | "Mol-weight of shoot reserve"
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

###########################################################################################
# Variables


" Model variables "
@plottable @units @udefault_kw struct Vars{F,MoMoD,C,M,T}
    shape::F             | [0.0]   | _                 | _
    rate::MoMoD          | [0.0]   | mol*mol^-1*d^-1   | _
    θE::F                | [0.0]   | _                 | _
    temp::C              | [25.0]  | K                 | _
    tempcorrection::F    | [1.0]   | _                 | _
    height::M            | [0.0]   | m                 | _
    t::T                 | [1]     | _                 | false
end

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

" Basic model components. For a plants, organs might be roots, stem and leaves "
struct Organ{P,S,V,F,F1} <: AbstractOrgan
    params::P
    shared::S
    vars::V
    J::F
    J1::F1
end

" update vars and flux to record for time t "
define_organs(o, t) = define_organs(o.params, o.shared, o.records, t)
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
    J::F    | false
    J1::F1  | false 
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

"An organism, made up of organs"
@flattenable mutable struct Plant{P,S,R,E,ES,D} <: AbstractOrganism
    params::P             | true
    shared::S             | true
    records::R            | false
    environment::E        | false
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
        environment_start = setindex!(Array{typeof(1hr),0}(undef), 1hr),
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
