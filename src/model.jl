"""
Run a DEB organism model from DiffEq.
"""
(organism::Plant)(du, u, t::Number, organs) = begin
    dead(organism) && return

    # Make sure the parameters don't any physical laws
    check_params(organs)

    # Split state into separate organs
    ux = split_state(organs, u)

    # Set up variables for this timestep and the current state
    apply(update_tstep!, organs, t)
    apply(zero_flux!, organs)
    apply(update_height!, organs, ux)
    apply(update_shape!, organs, ux)
    apply_environment!(organism, organs, ux, t)

    # Run the model, tag the organism as dead if it breaks.
    if !debmodel!(organs, ux, environment(organism)) 
        du .= zero(eltype(du)) 
        set_dead!(organism, true)
        return
    end

    # Sum the flux matrix to the state change vector
    sum_flux!(du, organs)
    nothing
end

"""
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation methods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
debmodel!(organs::Tuple, u::Tuple, env) = begin
    # Quit if it dies 
    false in metabolism!(organs, u) && return false
    translocation!(organs)
    assimilation!(organs, u)
    true
end

"""
Metabolism is an identical process for all organs, with potentially
different parameters or area and rate functions.
"""
metabolism!(organs::Tuple, u) = apply(metabolism!, organs, u)
metabolism!(o::AbstractOrgan, u) = begin
    # Quit if it dies 
    catabolism!(o, u) || return false
    growth!(o, u)
    maturity!(o, u)
    maintenence!(o, u)
    feedback!(o, u)
    true
end

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not alter flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = catabolism!(catabolism_pars(o), o, u)
catabolism!(p::CatabolismCNE, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (p.k_EC, p.k_EN, p.k_E) .* tempcorrection(v) .* shape(v)
    reserve = (u.C, u.N, u.E)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    set_var!(o, :rate, r)
    r < zero(r) && return false

    J1[:C,:ctb], J1[:N,:ctb], J1[:EE,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:CN,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))

    J1[:E,:ctb] = J1[:EE,:ctb] + J1[:CN,:ctb] # Total catabolic flux
    set_var!(o, :θE, J1[:EE,:ctb]/J1[:E,:ctb]) # Proportion of general reserve flux in total catabolic flux
    true
end
catabolism!(p::CatabolismCN, o, u) = begin
    v, J, J1 = vars(o), flux(o), flux1(o)
    turnover = (p.k_EC, p.k_EN) .* tempcorrection(v) .* shape(v)
    reserve = (u.C, u.N)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    set_var!(o, :rate, r)
    r < zero(r) && return false

    J1[:C,:ctb], J1[:N,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:E,:ctb] = stoich_merge(su_pars(o), J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    true
end

"""
Allocates reserves to growth.
"""
growth!(o, u) = begin
    flux(o)[:V,:gro] = growth = rate(o) * u.V
    drain = (1/y_V_E(o)) * growth 
    product = growth_production!(o, growth)
    reserve_drain!(o, Val(:gro), drain)
    # loss = drain - growth - product
    # reserve_loss!(o, loss)
    # conversion_loss!(o, growth, n_N_V(o))
end

growth_production!(o, growth) = growth_production!(production_pars(o), o, growth)
growth_production!(p::Production, o, growth) = flux(o)[:P,:gro] = growth * p.y_P_V
growth_production!(p, o, growth) = zero(growth)

"""
Allocates reserve drain due to maintenance.
"""
maintenence!(o, u) = begin
    drain = j_E_mai(o) * tempcorrection(o) * u.V
    maint_prod = maintenance_production!(o, u)
    reserve_drain!(o, Val(:mai), drain)
    # reserve_loss!(o, drain - maint_prod) # all maintenance is loss
end

maintenance_production!(o, u) = maintenance_production!(production_pars(o), o, u)
maintenance_production!(p::Production, o, u) = flux(o)[:P,:mai] = p.j_P_mai * tempcorrection(o) * u.V
maintenance_production!(p, o, u) = zero(eltype(flux(o)))

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Maturity, o, u) = begin
    flux(o)[:M,:gro] = mat = f.κmat * flux1(o)[:E,:ctb]
    mat_mai = f.j_E_mat_mai * tempcorrection(v) * u.V # min(u[:V], f.M_Vmat))
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
    # reserve_loss!(o, mat_mai)
    # conversion_loss!(o, mat, f.n_N_M)
end
maturity!(f::Nothing, o, u) = nothing

"""
    assimilation!(o, u)
Runs assimilation methods, depending on formulation and state.
"""
assimilation!(organs::Tuple, u) = apply(assimilation!, organs, u)
assimilation!(o::AbstractOrgan, u) = begin
    is_germinated(o, u) && assimilation!(has_reserves(o), assimilation_pars(o), o, u)
    nothing
end
assimilation!(::Nothing, x, o::AbstractOrgan, u, env) = nothing

"""
Translocation occurs between adjacent organs.
This function is identical both directiono, so on represents
whichever is not the current organs.

Will not run with less than 2 organs.
"""
translocation!(organs::Tuple{Organ, Organ}) = begin
    translocate_rejected!(organs[1], organs[2], 1.0)
    translocate_rejected!(organs[2], organs[1], 1.0)
    translocate!(organs[1], organs[2], 1.0)
    translocate!(organs[2], organs[1], 1.0)
end
translocation!(organs::Tuple) = translocation!(organs...)
translocation!(organs::Tuple{}) = nothing
translocation!(organs::Tuple{Organ}) = nothing


"""
Reallocate state rejected from synthesizing units.
TODO: add a 1-organs method. How does this interact with assimilation?  
"""
translocate_rejected!(source, dest, prop) = translocate_rejected!(rejection_pars(source), source, dest, prop)
translocate_rejected!(rejected::Nothing, source, dest, prop) = nothing
translocate_rejected!(rejected::DissipativeRejection, source, dest, prop) = begin
    Js, J1s, Jd = flux(source), flux1(source), flux(dest)
    transC = J1s[:C,:rej] # * (1 - κEC(o))
    transN = J1s[:N,:rej] # * (1 - κEN(o))
    Js[:C,:rej] = -transC
    Js[:N,:rej] = -transN
    Jd[:C,:tra] = y_EC_ECT(o) * transC
    Jd[:N,:tra] = y_EN_ENT(o) * transN
    # J1[:C,:los] += transC * (1 - y_EC_ECT(o)) + transN * (1 - y_EN_ENT(o))
    # J1[:N,:los] += (transC * (1 - y_EC_ECT(o)), transN * (1 - y_EN_ENT(o))) ⋅ (n_N_EC(o), n_N_EN(o))
    nothing
end
translocate_rejected!(rejected::LosslessRejection, source, dest, prop) = begin
    Js, J1s, Jd = flux(source), flux1(source), flux(dest)
    transC = J1s[:C,:rej] # * (1 - κEC(o))
    transN = J1s[:N,:rej] # * (1 - κEN(o))
    Js[:C,:rej] = -transC
    Js[:N,:rej] = -transN
    Jd[:C,:tra] = transC
    Jd[:N,:tra] = transN
    nothing
end

"""
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
@inline reserve_drain!(o::AbstractOrgan, col, drain) = reserve_drain!(has_reserves(o), o, col, drain)
@inline reserve_drain!(::HasCNE, o, col, drain) = begin
    ΘE, J = θE(o), flux(o)
    J_CN = -drain * (1 - θ) # fraction on drain from C and N reserves
    @inbounds J[:C,col] = J_CN/y_E_EC(o)
    @inbounds J[:N,col] = J_CN/y_E_EN(o)
    @inbounds J[:E,col] = -drain * ΘE 
    nothing
end
@inline reserve_drain!(::HasCN, o, col, drain) = begin
    J = flux(o)
    @inbounds J[:C,col] = -drain/y_E_EC(o)
    @inbounds J[:N,col] = -drain/y_E_EN(o)
    nothing
end

update_shape!(o, u) = set_shape!(o, shape_correction(shape_pars(o), u.V))

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
is_germinated(o, u) = is_germinated(germination_pars(o), o, u)
is_germinated(g::Nothing, o, u) = true
is_germinated(g::Germination, o, u) = u.V > g.M_Vgerm


update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(a::Nothing, o, u) = nothing
update_height!(a, o, u) = set_height!(o, allometric_height(a, o, u))

allometric_height(f::SqrtAllometry, o, u) = begin
    units = oneunit(u.V * w_V(o))
    sqrt((u.V * w_V(o)) / units) * f.allometry
end

# J: Flux matrix diagram.
# Rows: state.
# Columns: transformations
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ assS       │ groS       │ maiS       │ matS       │ rejS       │ traS       ┃ assR       │ groR       │ maiR       │ ratR │ rejR       │ traR       ┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃    ┃ JSS SubArray                                                                ┃ JRR SubArray                                                          ┃
# ┃    ┃                                                                             ┃                                                                       ┃
# ┃PS  ┃ 0          │ J_PS,groS  │ J_PS,maiS  │ 0          │ 0          │ 0          ┃ 0          │ J_PR,groR  │ J_PR,maiR  │ 0    │ 0          │ 0          ┃
# ┃VS  ┃ 0          │ J_VS,groS  │ 0          │ 0          │ 0          │ 0          ┃ 0          │ J_VR,groR  │ 0          │ 0    │ 0          │ 0          ┃
# ┃RS  ┃ 0          │ 0          │ 0          │ J_MS,groS  │ 0          │ 0          ┃ 0          │ 0          │ 0          │ 0    │ 0          │ 0          ┃
# ┃ECS ┃ J_ECS,assS │ J_ECS,groS │ J_ECS,maiS │ J_ECS,matS │ J_ECS,rejS │ J_ECS,traS ┃ J_ECR,assR │ J_ECR,groR │ J_ECR,maiR │ 0    │ J_ECS,rejR │ J_ECR,traR ┃
# ┃ENS ┃ J_ENS,assS │ J_ENS,groS │ J_ENS,maiS │ J_ENS,matS │ J_ENS,rejS │ J_ENS,traS ┃ J_ENR,assR │ J_ENR,groR │ J_ENR,maiR │ 0    │ J_ENS,rejR │ J_ENR,traR ┃
# ┃ES  ┃ J_ES,assS  │ J_ES,groS  │ J_ES,maiS  │ J_ES,matS  │ 0          │ J_ES,traS  ┃ J_ER,assR  │ J_ER,groR  │ J_ER,maiR  │ 0    │ 0          │ J_ER,traR  ┃
# ┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# J1: Catabolic flux diagram.
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ catS       │ rejS      │ losS mols!┃    ┃ catR       │ rejR      │ losR mols!┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃EES ┃ J_EES,catS │ 0         │ 0         ┃EES ┃ J_EER,catR │ 0         │ 0         ┃
# ┃CNS ┃ J_CNS,catS │ 0         │ 0         ┃NS  ┃ J_CNR,catR │ 0         │ 0         ┃
# ┃CS  ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃CS  ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
# ┃NS  ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃ENS ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
# ┃ES  ┃ J_ES,catS  │ 0         │ 0         ┃ES  ┃ J_ER,catR  │ 0         │ 0         ┃
# ┗━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# M: State vector diagram.
# ┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
