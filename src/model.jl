"""
Run a DEB organism model.
"""
(o::Organism)(du, u, p::Nothing, t::Number) = begin
    update_t!(o.organs, t)
    o(du, u, t, o.organs)
end

(o::Organism)(du, u, t::Number, organs) = begin
    check_params(organs)
    o.dead[] && return
    ux = split_state(organs, u)
    apply(set_height!, organs, ux)
    apply_environment!(organs, ux, o.environment, t)
    if !debmodel!(organs, ux) 
        du .= zero(eltype(du)) 
        o.dead[] = true
        return
    end
    sum_flux!(du, organs)
    return
end

"""
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
debmodel!(organs::Tuple, u::Tuple) = begin
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
metabolism!(o::Organ, u) = begin
    set_scale!(o, u)
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
Does not finalise flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = catabolism!(turnover_pars(o), o, u)
catabolism!(t::TurnoverCNE, o, u) = begin
    v, J, J1 = unpack(o)
    turnover = (t.k_EC, t.k_EN, t.k_E) .* tempcorrection(v) .* scale(v)
    reserve = (u.C, u.N, u.E)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    set_var!(o, :rate, r)
    r < zero(r) && return false

    J1[:C,:ctb], J1[:N,:ctb], J1[:EE,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:CN,:ctb] = stoich_merge(o, J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))

    J1[:E,:ctb] = J1[:EE,:ctb] + J1[:CN,:ctb] # Total catabolic flux
    set_var!(o, :θE, J1[:EE,:ctb]/J1[:E,:ctb]) # Proportion of general reserve flux in total catabolic flux
    true
end
catabolism!(t::TurnoverCN, o, u) = begin
    v, J, J1 = unpack(o)
    turnover = (t.k_EC, t.k_EN) .* tempcorrection(v) .* scale(v)
    reserve = (u.C, u.N)
    rel_reserve = reserve ./ u.V
    corr_j_E_mai = j_E_mai(o) * tempcorrection(v)

    r = calc_rate(rate_formula(o), su_pars(o), rel_reserve, turnover, corr_j_E_mai, y_E_EC(o), y_E_EN(o), y_V_E(o), κsoma(o))
    set_var!(o, :rate, r)
    r < zero(r) && return false

    J1[:C,:ctb], J1[:N,:ctb] = non_growth_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:E,:ctb] = stoich_merge(o, J1[:C,:ctb], J1[:N,:ctb], y_E_EC(o), y_E_EN(o))
    true
end

"""
Allocates reserves to growth.
"""
growth!(o, u) = begin
    drain = rate(o.vars) * u.V
    o.J[:V,:gro] = growth = y_V_E(o) * drain
    product = growth_production!(o, growth)
    reserve_drain!(o, Val(:gro), drain)
    # loss = drain - growth - product
    # reserve_loss!(o, loss)
    # conversion_loss!(o, growth, n_N_V(o))
end

growth_production!(o, growth) = growth_production!(production_pars(o), o, growth)
growth_production!(p::Production, o, growth) = o.J[:P,:gro] = growth * p.y_P_V
growth_production!(p, o, growth) = zero(growth)

"""
Allocates reserve drain due to maintenance.
"""
maintenence!(o, u) = begin
    drain = j_E_mai(o) * tempcorrection(o.vars) * u.V
    maint_prod = maintenance_production!(o, u)
    reserve_drain!(o, Val(:mai), drain)
    # reserve_loss!(o, drain - maint_prod) # all maintenance is loss
end

maintenance_production!(o, u) = maintenance_production!(production_pars(o), o, u)
maintenance_production!(p::Production, o, u) = o.J[:P,:mai] = p.j_P_mai * tempcorrection(o.vars) * u.V
maintenance_production!(p, o, u) = zero(eltype(o.J))

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Maturity, o, u) = begin
    o.J[:M,:gro] = mat = f.κmat * o.J1[:E,:ctb]
    mat_mai = f.j_E_mat_mai * tempcorrection(v) * u.V # min(u[:V], f.M_Vmat))
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
    # reserve_loss!(o, mat_mai)
    # conversion_loss!(o, mat, f.n_N_M)
end
maturity!(f::Nothing, o, u) = nothing

"""
Translocation occurs between adjacent organs.
This function is identical both directiono, so on represents
whichever is not the current organs.

Will not run with less than 2 organs.
"""
translocation!(organs::Tuple{Organ, Organ}) = begin
    reuse_rejected!(organs[1], organs[2], 1.0)
    reuse_rejected!(organs[2], organs[1], 1.0)
    translocate!(organs[1], organs[2], 1.0)
    translocate!(organs[2], organs[1], 1.0)
end
translocation!(organs::Tuple) = translocation!(organs...)
translocation!(organs::Tuple{}) = nothing
translocation!(organs::Tuple{Organ}) = nothing

# Recurse through all organs. A loop would not be type-stable.
# translocation!(organs::Tuple, destorgans::Tuple) = begin
#     props = buildprops(organs[1])
#     translocation!(organs[1], destorgans, organs[1].params.translocation.destnames, props)
#     translocation!(tail(organs), destorgans)
# end
# translocation!(organs::Tuple{}, destorgans::Tuple) = nothing
# translocation!(organ::Organ, destorgans::Tuple, destnames::Symbol, props) =
#     translocation!(organ, destorgans, (destnames,), props)
# # Translocate to organs with names in the destnames list
# translocation!(organ::Organ, destorgans::Tuple, destnames, props) = begin
#     for i = 1:length(destnames)
#         if destorgans[1].params.name == destnames[i]
#             reuse_rejected!(organ, destorgans[1], props[i])
#             translocate!(organ, destorgans[1], props[i])
#             break
#         end
#     end
#     translocation!(organ, tail(destorgans), destnames, props)
# end
# translocation!(organ::Organ, destorgans::Tuple{}, destnames, props) = nothing

# Add the last remainder proportion (so that its not a model parameter)
# buildprops(o::Organ) = buildprops(o.params.translocation.proportions)
# buildprops(x::Nothing) = (1.0)
# buildprops(x::Number) = (x, 1 - x)
# buildprops(xs::Tuple) = (xs..., 1 - sum(xs))

"""
Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs.
This function is identical both directiono, and ox represents
whichever is not the current organs. Will not run with less than 2 organs.
"""
translocate!(o1, o2, prop) = translocate!(trans_pars(o1), o1, o2, prop)
translocate!(p::Nothing, o1, o2, prop) = nothing
translocate!(p::AbstractDissipativeTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx * y_E_ET(o2)

    # loss = transx * (1 - y_E_ET(o2))
    # reserve_loss!(o2, loss)
    # conversion_loss!(o2, transx * y_E_ET(o2), n_N_E(o2))
    nothing
end

translocate!(p::AbstractLosslessTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx

    # conversion_loss!(o2, transx, n_N_E(o2))
    nothing
end

"""
Reallocate state rejected from synthesizing units.
TODO add a 1-organs method Also how does this interact with assimilation?  """
reuse_rejected!(source, dest, prop) = reuse_rejected!(rejection_pars(source), source, dest, prop)
reuse_rejected!(rejected::Nothing, source, dest, prop) = nothing
reuse_rejected!(rejected::DissipativeRejection, source, dest, prop) = begin
    v, J, J1 = unpack(o)
    transC = J1[:C,:rej] # * (1 - κEC(o))
    transN = J1[:N,:rej] # * (1 - κEN(o))
    J[:C,:rej] = -transC
    J[:N,:rej] = -transN
    dest.J[:C,:tra] = y_EC_ECT(o) * transC
    dest.J[:N,:tra] = y_EN_ENT(o) * transN
    J1[:C,:los] += transC * (1 - y_EC_ECT(o)) + transN * (1 - y_EN_ENT(o))
    J1[:N,:los] += (transC * (1 - y_EC_ECT(o)), transN * (1 - y_EN_ENT(o))) ⋅ (n_N_EC(o), n_N_EN(o))
    nothing
end
reuse_rejected!(rejected::LosslessRejection, source, dest, prop) = begin
    transC = source.J1[:C,:rej] # * (1 - κEC(o))
    transN = source.J1[:N,:rej] # * (1 - κEN(o))
    source.J[:C,:rej] = -transC
    source.J[:N,:rej] = -transN
    dest.J[:C,:tra] = transC
    dest.J[:N,:tra] = transN
    nothing
end

"""
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
@inline reserve_drain!(o::Organ, col, drain) = reserve_drain!(has_reserves(o), o, col, drain)
@inline reserve_drain!(::HasCNE, o, col, drain) = begin
    θ = θE(o)
    J_CN = -drain * (1 - θ) # fraction on drain from C and N reserves
    @inbounds o.J[:C,col] = J_CN/y_E_EC(o)
    @inbounds o.J[:N,col] = J_CN/y_E_EN(o)
    @inbounds o.J[:E,col] = -drain * θ
    nothing
end
@inline reserve_drain!(::HasCN, o, col, drain) = begin
    @inbounds o.J[:C,col] = -drain/y_E_EC(o)
    @inbounds o.J[:N,col] = -drain/y_E_EN(o)
    nothing
end

"""
Generalised reserve loss to track carbon.

Loss is distributed between general and C and N reserves by the fraction θE
"""
@inline reserve_loss!(o, loss) = nothing #reserve_loss!(o.params, o, loss)
@inline reserve_loss!(::HasCNE, o, loss) = begin
    ee = loss * θ # fraction of loss from E reserve
    ecn = loss - ee # fraction on loss from C and N reserves
    ec = ecn/y_E_EC(o)
    en = ecn/y_E_EN(o)
    o.J1[:C,:los] += ec + en + ee
    o.J1[:N,:los] += (ec, en, ee) ⋅ (n_N_EC(o), n_N_EN(o), n_N_E(o))
    nothing
end
@inline reserve_loss!(::HasCN, o, loss) = begin
    ec = loss/y_E_EC(o)
    en = loss/y_E_EN(o)
    o.J1[:C,:los] += ec + en
    o.J1[:N,:los] += (ec, en) ⋅ (n_N_EC(o), n_N_EN(o))
    nothing
end

@inline conversion_loss!(o, loss, dest_n_N) = nothing #conversion_loss!(o.params, o, loss, dest_n_N)
@inline conversion_loss!(::HasCNE, o, loss, dest_n_N) = begin
    ecn = loss * (1 - θE(o)) # fraction on loss from C and N reserves
    ec = ecn/y_E_EC(o)
    en = ecn/y_E_EN(o)
    o.J1[:C,:los] += ec + loss * (θE(o) - 1) # + en
    o.J1[:N,:los] += (ec, en, loss * (θ - dest_n_N/n_N_E(o))) ⋅ (n_N_EC(o), n_N_EN(o), n_N_E(o))
end
@inline conversion_loss!(::HasCN, o, loss, dest_n_N) = begin
    ec = loss/y_E_EC(o)
    en = loss/y_E_EN(o)
    o.J1[:C,:los] += ec + loss # + en
    o.J1[:N,:los] += (ec, en) ⋅ (n_N_EC(o), n_N_EN(o))
end

"""
    stoich_merge(Ja, Jb, ya, yb)
Merge fluxes stoichiometrically into general reserve Eab based on yeild
fractions ya and yb. An unmixed proportion is returned as unmixed reserves Ea and Eb.
Losses are also calculated.
"""
stoich_merge(o, Ja, Jb, ya, yb) = begin
    JEab = synthesizing_unit(su_pars(o), Ja * ya, Jb * yb)
    Ja1 = Ja - JEab/ya
    Jb1 = Jb - JEab/yb
    (Ja1, Jb1, JEab)
end

stoich_merge_losses(Jc1, Jn1, Jc2, Jn2, JEcn, n_c, n_n, n_Ecn) = begin
    lossa = Jc1 - Jc2 + Jn1 - Jn2 - JEcn
    lossb = (Jc1 - Jc2, Jn1 - Jn2, -JEcn) ⋅ (n_c, n_n, n_Ecn)
    lossa, lossb
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.
"""
feedback!(o, u) = feedback!(feedback_pars(o), has_reserves(o), o, u)
feedback!(f::Nothing, x, o, u) = nothing
feedback!(f::Autophagy, ::HasCNE, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    # TODO this should be a lossy process
    o.J[:E,:fbk] += aph
    o.J[:V,:fbk] -= aph
    nothing
end
feedback!(f::Autophagy, ::HasCN, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    o.J[:C,:fbk] += aph
    o.J[:N,:fbk] += aph * n_N_V(o)
    o.J[:V,:fbk] -= aph
    nothing
end


set_height!(o, u) = set_height!(allometry_pars(o), o, u)
set_height!(a::Nothing, o, u) = nothing
set_height!(a, o, u) = begin
    h = allometric_height(a, o, u)
    set_var!(o, :height, h)
end

set_scale!(o, u) = set_var!(o, :scale, calc_scaling(scaling_pars(o), u.V))

calc_scaling(f::KooijmanArea, uV) = begin
    # uV <= zero(uV) && error("Mass is less than zero, I think its dead...")
    (uV / f.M_Vref)^(-uV / f.M_Vscaling)
end
calc_scaling(f::Nothing, uV) = 1

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
is_germinated(o, u) = u.V > M_Vgerm(o)


allometric_height(o, u) = allometric_height(allometry_pars(o), o, u)
allometric_height(f::Nothing, o, u) = zero(height(o.vars))
allometric_height(f::SqrtAllometry, o, u) = begin
    dim = oneunit(u.V * w_V(o))
    sqrt((u.P * w_P(o) + u.V * w_V(o)) / dim) * f.allometry
end

unpack(o::Organ) = o.vars, o.J, o.J1

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
