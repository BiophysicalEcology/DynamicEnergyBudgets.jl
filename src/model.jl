"""
Run a DEB organism model.
"""
(o::Organism)(du, u, p::Nothing, t::Number) = o(du, u, t, define_organs(o, t))

(o::Organism)(du, u, t::Number, organs) = begin
    check_params.(organs)
    ux = split_state(organs, u)
    set_height!.(organs, ux)
    # apply_environment!(organs, ux, o.environment, t)
    debmodel!(organs, ux)
    sum_flux!(du, organs)
end

set_height!(o, u) = set_height!(o.params.allometry, o, u)
set_height!(a::Nothing, o, u) = nothing
set_height!(a, o, u) = begin
    h = allometric_height(a, o, u)
    set_var!(o.vars, :height, h)
end

"""
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
debmodel!(organs::Tuple, u::Tuple) = begin
    metabolism!(organs, u)
    translocation!(organs)
    assimilation!(organs, u)
end

"""
Metabolism is an identical process for all organs, with potentially
different parameters or area and rate functions.
"""
metabolism!(organs::Tuple, u) = apply(metabolism!, organs, u)
metabolism!(o, u) = begin
    set_scaling!(o, u)
    catabolism!(o, u)
    growth!(o, u)
    maturity!(o, u)
    maintenence!(o, u)
    feedback!(o, u)
end

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not finalise flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = catabolism!(o.params, o, u)
catabolism!(p::ParamsCNE, o, u) = begin
    _, v, sh, J, J1 = unpack(o)
    turnover = (p.k_EC, p.k_EN, p.k_E) .* tempcorrection(v) .* scale(v)
    reserve = (u.C, u.N, u.E)
    rel_reserve = reserve ./ u.V
    j_E_mai = p.j_E_mai * tempcorrection(v)

    r = find_rate(rel_reserve, turnover, j_E_mai, p.y_E_CH_NO, p.y_E_EN, p.y_V_E, κsoma(o))
    set_var!(v, :rate, r)

    J1[:C,:ctb], J1[:N,:ctb], J1[:EE,:ctb] = catabolic_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:CN,:ctb] = stoich_merge(J1[:C,:ctb], J1[:N,:ctb], p.y_E_CH_NO, p.y_E_EN)

    J1[:E,:ctb] = J1[:EE,:ctb] + J1[:CN,:ctb] # Total catabolic flux
    set_var!(v, :θE, J1[:EE,:ctb]/J1[:E,:ctb]) # Proportion of general reserve flux in total catabolic flux
    nothing
end
catabolism!(p::ParamsCN, o, u) = begin
    _, v, sh, J, J1 = unpack(o)
    turnover = (p.k_EC, p.k_EN) .* tempcorrection(v) .* scale(v)
    reserve = (u.C, u.N)
    rel_reserve = reserve ./ u.V
    j_E_mai = p.j_E_mai * tempcorrection(v)

    r = find_rate(rel_reserve, turnover, j_E_mai, p.y_E_CH_NO, p.y_E_EN, p.y_V_E, κsoma(o))
    set_var!(v, :rate, r)

    J1[:C,:ctb], J1[:N,:ctb] = catabolic_flux.(reserve, turnover, r)
    J1[:C,:rej], J1[:N,:rej], J1[:E,:ctb] = stoich_merge(J1[:C,:ctb], J1[:N,:ctb], p.y_E_CH_NO, p.y_E_EN)
    nothing
end


"""
Allocates reserves to growth.
"""
growth!(o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    J[:V,:gro] = growth = rate(v) * u.V
    J[:P,:gro] = production = growth * p.y_P_V
    drain = (1/p.y_V_E) * growth 
    loss = drain - growth - production
    reserve_drain!(o, Val(:gro), drain)
    reserve_loss!(o, loss)
    conversion_loss!(o, growth, sh.n_N_V)
end

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(o.params.maturity, o, u)
maturity!(f::Maturity, o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    J[:M,:gro] = maturity = f.κmat * J1[:E,:ctb]
    mat_mai = f.j_E_mat_mai * tempcorrection(v) * u.V # min(u[:V], f.M_Vmat))
    drain = maturity + mat_mai 
    reserve_drain!(o, Val(:mat), drain)
    reserve_loss!(o, mat_mai)
    conversion_loss!(o, maturity, f.n_N_M)
end
maturity!(f::Nothing, o, u) = nothing

"""
Allocates reserve drain due to maintenance.
"""
maintenence!(o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    # TODO this isn't balanced anywhere and is
    # larger than the total maintenance flux
    drain = p.j_E_mai * tempcorrection(v) * u.V
    J[:P,:mai] = maint_prod = p.j_P_mai * tempcorrection(v) * u.V
    reserve_drain!(o, Val(:mai), drain)
    reserve_loss!(o, drain - maint_prod) # all maintenance is loss
end

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
translocate!(o1, o2, prop) = translocate!(o1.params.translocation, o1, o2, prop)
translocate!(p::Nothing, o1, o2, prop) = nothing
translocate!(p::AbstractDissipativeTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx * o2.params.y_E_ET

    loss = transx * (1 - o2.params.y_E_ET)
    reserve_loss!(o2, loss)
    conversion_loss!(o2, transx * o2.params.y_E_ET, o2.shared.n_N_E)
    nothing
end

translocate!(p::AbstractLosslessTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx

    conversion_loss!(o2, transx, o2.shared.n_N_E)
    nothing
end

"""
Reallocate state rejected from synthesizing units.
TODO add a 1-organs method Also how does this interact with assimilation?  """
reuse_rejected!(source, dest, prop) = reuse_rejected!(source.params.rejection, source, dest, prop)
reuse_rejected!(rejected::Nothing, source, dest, prop) = nothing
reuse_rejected!(rejected::DissipativeRejection, source, dest, prop) = begin
    p, v, sh, J, J1 = unpack(source);

    transC = J1[:C,:rej] # * (1 - p.κEC)
    transN = J1[:N,:rej] # * (1 - p.κEN)
    J[:C,:rej] = -transC 
    J[:N,:rej] = -transN
    dest.J[:C,:tra] = p.y_EC_ECT * transC
    dest.J[:N,:tra] = p.y_EN_ENT * transN
    J1[:C,:los] += transC * (1 - p.y_EC_ECT) + transN * (1 - p.y_EN_ENT)
    J1[:N,:los] += (transC * (1 - p.y_EC_ECT), transN * (1 - p.y_EN_ENT)) ⋅ (sh.n_N_EC, sh.n_N_EN)
    nothing
end
reuse_rejected!(rejected::LosslessRejection, source, dest, prop) = begin
    p, v, sh, J, J1 = unpack(source);

    transC = J1[:C,:rej] # * (1 - p.κEC)
    transN = J1[:N,:rej] # * (1 - p.κEN)
    J[:C,:rej] = -transC 
    J[:N,:rej] = -transN
    dest.J[:C,:tra] = transC
    dest.J[:N,:tra] = transN
    nothing
end

"""
Generalised reserve drain for any flux column *col* (ie :gro)
and any combination of reserves.
"""
@inline reserve_drain!(o::Organ, col, drain) = reserve_drain!(o.params, o, col, drain)
@inline reserve_drain!(p::ParamsCNE, o, col, drain) = begin
    θ = θE(o)
    J_CN = -drain * (1 - θ) # fraction on drain from C and N reserves
    @inbounds o.J[:C,col] = J_CN/o.params.y_E_CH_NO
    @inbounds o.J[:N,col] = J_CN/o.params.y_E_EN
    @inbounds o.J[:E,col] = -drain * θ
    nothing
end
@inline reserve_drain!(p::ParamsCN, o, col, drain) = begin
    @inbounds o.J[:C,col] = -drain/o.params.y_E_CH_NO
    @inbounds o.J[:N,col] = -drain/o.params.y_E_EN
    nothing
end

"""
Generalised reserve loss to track carbon. 

Loss is distributed between general and C and N reserves by the fraction θE
"""
reserve_loss!(o, loss) = nothing #reserve_loss!(o.params, o, loss)
reserve_loss!(p::ParamsCNE, o, loss) = begin
    p, v, sh = unpack(o); θ = θE(o)
    ee = loss * θ # fraction of loss from E reserve
    ecn = loss - ee # fraction on loss from C and N reserves
    ec = ecn/p.y_E_CH_NO
    en = ecn/p.y_E_EN
    o.J1[:C,:los] += ec + en + ee
    o.J1[:N,:los] += (ec, en, ee) ⋅ (sh.n_N_EC, sh.n_N_EN, sh.n_N_E)
    nothing
end
reserve_loss!(p::ParamsCN, o, loss) = begin
    p, v, sh = unpack(o)
    ec = loss/p.y_E_CH_NO
    en = loss/p.y_E_EN
    o.J1[:C,:los] += ec + en
    o.J1[:N,:los] += (ec, en) ⋅ (sh.n_N_EC, sh.n_N_EN)
    nothing
end

conversion_loss!(o, loss, dest_n_N) = nothing #conversion_loss!(o.params, o, loss, dest_n_N)
conversion_loss!(p::ParamsCNE, o, loss, dest_n_N) = begin
    p, v, sh = unpack(o); θ = θE(o)
    ecn = loss * (1 - θ) # fraction on loss from C and N reserves
    ec = ecn/p.y_E_CH_NO
    en = ecn/p.y_E_EN
    o.J1[:C,:los] += ec + loss * (θ - 1) # + en
    o.J1[:N,:los] += (ec, en, loss * (θ - dest_n_N/sh.n_N_E)) ⋅ (sh.n_N_EC, sh.n_N_EN, sh.n_N_E)
end
conversion_loss!(p::ParamsCN, o, loss, dest_n_N) = begin
    p, v, sh = unpack(o)
    ec = loss/p.y_E_CH_NO
    en = loss/p.y_E_EN
    o.J1[:C,:los] += ec + loss # + en
    o.J1[:N,:los] += (ec, en) ⋅ (sh.n_N_EC, sh.n_N_EN)
end

"""
    stoich_merge(Ja, Jb, ya, yb)
Merge fluxes stoichiometrically into general reserve Eab based on yeild
fractions ya and yb. An unmixed proportion is returned as unmixed reserves Ea and Eb.
Losses are also calculated.
"""
stoich_merge(Ja, Jb, ya, yb) = begin
    JEab = synthesizing_unit(Ja * ya, Jb * yb) 
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
κtra is the difference paramsbetween κsoma and κmat
"""
κtra(o::Organ) = κtra(o.params)
κtra(p) = κtra(p.translocation)
κtra(translocation::AbstractTranslocation) = translocation.κtra
κtra(o::Nothing) = 0.0

κmat(o::Organ) = κmat(o.params)
κmat(p) = κmat(p.maturity)
κmat(maturity::Maturity) = maturity.κmat
κmat(maturity::Nothing) = 0.0

κsoma(o::Organ) = (1.0 - κtra(o) - κmat(o))

"""
Calculate rate formula. TODO: use Roots.jl for this
"""
function find_rate(m, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)
    # bounds = (-10.0oneunit(v.rate), 20.0oneunit(v.rate)) #rate_window(args...)

    # Find the type of rate so that diff and units work with the guess and atol
    one_r = oneunit(eltype(turnover))
    # Get rate with a zero finder
    let m=m, turnover=turnover, j_E_mai=j_E_mai, y_E_CH_NO=y_E_CH_NO, y_E_EN=y_E_EN, y_V_E=y_V_E, κsoma=κsoma
        find_zero((-2one_r, 1one_r), Secant(), one_r*1e-10, 100) do x
            rate_formula(x, m, turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)
        end
    end
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.
"""
feedback!(o, u) = feedback!(o.shared.feedback, o.params, o, u)

feedback!(f::Nothing, x, o, u) = nothing
feedback!(f::Autophagy, p::ParamsCNE, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    # TODO this should be a lossy process
    o.J[:E,:fbk] += aph
    o.J[:V,:fbk] -= aph
    nothing
end
feedback!(f::Autophagy, p::ParamsCN, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u.V * (oneunit(hs) - hs)
    o.J[:C,:fbk] += aph # TODO divide this by y_CH_NO etc
    o.J[:N,:fbk] += aph
    o.J[:V,:fbk] -= aph
    nothing
end


set_scaling!(o, u) = set_var!(o.vars, :scale, scaling(o.params.scaling, u.V))

scaling(f::KooijmanArea, uV) = begin
    uV > zero(uV) || return 1.0 # || error("Mass is less than zero, I think its dead...")
    (uV / f.M_Vref)^(-uV / f.M_Vscaling)
end
scaling(f::Nothing, uV) = 1.0

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
germinated(o, u) = u.V > o.params.M_Vgerm


allometric_height(f::SqrtAllometry, o, u) = begin
    p, v, sh = unpack(o)
    dim = oneunit(u.V * sh.w_V)
    sqrt((u.P * sh.w_P + u.V * sh.w_V) / dim) * f.allometry
end

unpack(o::Organ) = o.params, o.vars, o.shared, o.J, o.J1

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
