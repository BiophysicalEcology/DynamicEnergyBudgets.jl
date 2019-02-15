abstract type AbstractCatabolism end

@mix @columns struct CatabolismCN{MoMoD}
    k_EC::MoMoD | 0.4 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "C-reserve turnover rate"
    k_EN::MoMoD | 0.4 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "N-reserve turnover rate"
end

@mix @columns struct CatabolismE{MoMoD}
    k_E::MoMoD  | 0.4 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0,1.0]   | _  | "C-reserve turnover rate"
end

@CatabolismE struct CatabolismE{} <: AbstractCatabolism end
@CatabolismE struct CatabolismCN{} <: AbstractCatabolism end
@CatabolismE @CatabolismCN struct CatabolismCNE{} <: AbstractCatabolism end

# Unnecessary and probably suboptimal
# κEC::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]  | _ | "Non-processed C-reserve returned to C-reserve"
# κEN::F             | 0.3             | _               | Beta(2.0, 2.0)  | [0.0,1.0]  | _ | "Non-processed N-reserve returned to N-reserve"


abstract type AbstractMaintenance end

@columns struct Maintenance{MoMoD} <: AbstractMaintenance
    j_E_mai::MoMoD       | 0.023 | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [1e-4, 1.0] | true | "Spec somatic maint costs."
end


abstract type AbstractParams end

" Model parameters that vary between organs "
@default_kw @selectable struct Params{A,S,Al,Ma,Tr,Re,Ge,Pr} <: AbstractParams
    # Field              | Selectable Types                     | Default
    name::Symbol         | _                                    | :organ
    rate_formula         | _                                    | FZeroRate()
    assimilation_pars::A | Union{Nothing,AbstractAssimilation}  | ConstantCAssim()
    shape_pars::S        | Union{Nothing,AbstractShape}         | Plantmorph()
    allometry_pars::Al   | Union{Nothing,AbstractAllometry}     | nothing
    maturity_pars::Ma    | Union{Nothing,AbstractMaturity}      | nothing
    trans_pars::Tr       | Union{Nothing,AbstractTranslocation} | nothing
    rejection_pars::Re   | Union{Nothing,AbstractRejection}     | LosslessRejection()
    germination_pars::Ge | Union{Nothing,AbstractGermination}   | Germination()
    production_pars::Pr  | Union{Nothing,AbstractProduction}    | Production()
end


abstract type AbstractCore end

@columns struct Core{MoMo,GMo} <: AbstractCore
    y_V_E::MoMo  | 0.7   | mol*mol^-1 | Beta(2.0, 2.0)  | [0.0, 1.0]   | _ | "From reserve to structure"
    y_E_EC::MoMo | 0.7   | mol*mol^-1 | Gamma(2.0, 2.0) | [1e-6, 1.0]  | true | "From C-reserve to reserve, using nitrate"
    y_E_EN::MoMo | 30.0  | mol*mol^-1 | Gamma(2.0, 2.0) | [1.0, 50.0]  | false | "From N-reserve to reserve"
    n_N_V::MoMo  | 0.03  | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 0.1]   | _ | "N/C in structure"
    n_N_E::MoMo  | 0.025 | mol*mol^-1 | Gamma(2.0, 2.0) | [0.0, 0.1]   | _ | "N/C in reserve"
    w_V::GMo     | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [15.0, 40.0] | _ | "Mol-weight of shoot structure"
    w_N::GMo     | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [15.0, 40.0] | _ | "Mol-weight of shoot N-reserve"
    # w_C::GMo   | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [12.0, 40.0] | _ | "Mol-weight of shoot C-reserve"
    # w_E::GMo   | 25.0  | g*mol^-1   | Gamma(2.0, 2.0) | [15.0, 40.0] | _ | "Mol-weight of shoot reserve"
end

#    %   W   Rel n atoms
# C  45  12  30857
# H  6   1   58600
# O  45  16  27000
# N  2   14   1028
# K  1   39    246
# X  1   40    300ish


" Model parameters shared between organs "
@default_kw @selectable struct SharedParams{SU,Co,FB,TC,Tu,Mt}
    # Field              | Selectable Types                             | Default
    su_pars::SU          | Union{Nothing,AbstractSynthesizingUnit}      | ParallelComplementarySU()
    core_pars::Co        | nothing                                      | Core()
    feedback_pars::FB    | Union{Nothing,AbstractFeedback}              | nothing
    tempcorr_pars::TC    | Union{Nothing,AbstractTemperatureCorrection} | nothing
    catabolism_pars::Tu  | AbstractCatabolism                           | CatabolismCN()
    maintenance_pars::Mt | AbstractMaintenance                          | Maintenance()
end

###########################################################################################
# Variables

" Model variables "
@units @default_kw struct Vars{V,F,MoMoD,C,M,T}
    assimilation_vars::V | nothing | _
    shape::F             | [0.0]   | _
    rate::MoMoD          | [0.0]   | mol*mol^-1*d^-1
    θE::F                | [0.0]   | _
    temp::C              | [25.0]  | K
    tempcorrection::F    | [1.0]   | _
    height::M            | [0.0]   | m
    t::T                 | [1]     | _
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
struct Records{V,F,F1}
    vars::V
    J::F
    J1::F1
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
           vars = (ShootVars(), RootVars()),
           shared = SharedParams(),
           records = nothing,
           environment = nothing,
           time = 0hr:1hr:8760hr,
           environment_start = setindex!(Array{typeof(1hr),0}(undef), 1hr)
        ) = begin
    if records == nothing
        records = []
        for i = 1:length(params)
            push!(records, Records(params[i], vars[i], time, 1.0mol/hr, typeof(1.0mol/hr)))
        end
        records = (records...,)
    end
    dead = setindex!(Array{Bool,0}(undef), false)
    Plant(params, shared, records, environment, environment_start, dead)
end
