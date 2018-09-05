"""
Run a DEB organism model.
"""
(o::Organism)(du, u, p::Nothing, t::Number) = o(du, u, t, define_organs(o, t))

(o::Organism)(du, u, t::Number, organs) = begin
    ux = split_state(organs, u)
    apply_environment!(organs, ux, o.environment, t)
    debmodel!(organs, ux)
    sum_flux!(du, organs)
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
    feedback!(o.shared.feedback, o, u)
end

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not finalise flux in J - operates only on J1 (intermediate storage)
"""
catabolism!(o, u) = begin
    p, v, sh, J, J1 = unpack(o); va = assimilation(v);
    turnover = (p.k_EC, p.k_EN, p.k_E) .* tempcorrection(v) .* scale(v)
    reserve = (u[C], u[N], u[E])
    rel_reserve = reserve ./ u[V]
    j_E_mai = p.j_E_mai * tempcorrection(v)

    r = find_rate(rel_reserve, turnover, j_E_mai, sh.y_E_CH_NO, sh.y_E_EN, p.y_V_E, p.κsoma)
    set_var!(v, :rate, r)
    # r > zero(r) || error("rate is less than zero, thermodyamics disagrees")

    J1[C,ctb], J1[N,ctb], J1[EE,ctb] = catabolic_flux.(reserve, turnover, r)
    J1[C,rej], J1[N,rej], J1[CN,ctb] = stoich_merge(J1[C,ctb], J1[N,ctb], sh.y_E_CH_NO, sh.y_E_EN)

    J1[E,ctb] = J1[EE,ctb] + J1[CN,ctb] # Total catabolic flux
    set_var!(v, :θE, J1[EE,ctb]/J1[E,ctb]) # Proportion of general reserve flux in total catabolic flux
    nothing
end

"""
Allocates reserves to growth.
"""
growth!(o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    θ = θE(v)
    J[V,gro] = growth = rate(v) * u[V]
    drain = (1/p.y_V_E) * growth 
    loss = drain - growth
    reserve_drain!(o, gro, drain, θ)
    reserve_loss!(o, loss, θ)
    conversion_loss!(o, growth, θ, sh.n_N_V)
end

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(o.params.maturity, o, u)
maturity!(f::Maturity, o, u) = begin
    p, v, sh, J, J1 = unpack(o); θ = θE(v) 
    J[M,gro] = maturity = f.κmat * J1[E,ctb]
    mat_mai = f.j_E_mat_mai * tempcorrection(v) * u[V] # min(u[V], f.M_Vmat))
    drain = maturity + mat_mai 
    reserve_drain!(o, mat, drain, θ)
    reserve_loss!(o, mat_mai, θ)
    conversion_loss!(o, maturity, θ, sh.n_N_M)
end
maturity!(f::Nothing, o, u) = nothing

"""
Allocates reserve drain due to maintenance.
"""
maintenence!(o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    # J[P,mai] = maint_prod = p.j_P_mai * tempcorrection(v) * u[V]
    drain = p.j_E_mai * tempcorrection(v) * u[V]
    reserve_drain!(o, mai, drain, θE(v))
    reserve_loss!(o, drain, θE(v)) # all maintenance is loss
    # J1[C,los] -= maint_prod 
    # J1[N,los] -= maint_prod * sh.n_N_E
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
translocate!(o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[E,ctb]
    # reserve_drain!(o1, tra, trans, θE(o1.vars))

    # incoming translocation
    transx = κtra(o2) * o2.J1[E,ctb] * o2.params.y_E_ET 
    o1.J[E,tra] += transx * o2.params.y_E_ET
    # reserve_loss!(o2, transx * (1 - o2.params.y_E_ET), θE(o2.vars))
    conversion_loss!(o2, o1.J[E,tra], θE(o2.vars), o2.shared.n_N_V)
    nothing
end

"""
Reallocate state rejected from synthesizing units.

TODO add a 1-organs method
Also how does this interact with assimilation?
"""
reuse_rejected!(source, dest, prop) = begin
    p, v, sh, J, J1 = unpack(source);
    J[C,rej] = -J1[C,rej]
    J[N,rej] = -J1[N,rej]
    # Some rejected reserves are translocated and used in assimilation.
    # Why isn't the other rejected C and N used?
    if typeof(p.assimilation) <: AbstractCAssim
        dest.J[C,tra] = p.y_EC_ECT * J1[C,rej]
    elseif typeof(p.assimilation) <: AbstractNAssim
        dest.J[N,tra] = p.y_EN_ENT * J1[N,rej]
    end
    J1[C,los] = (1 - p.y_EC_ECT) * J1[C,rej] + (1 - p.y_EN_ENT) * J1[N,rej]
    J1[N,los] = (1 - p.y_EC_ECT) * J1[C,rej] * sh.n_N_EC + (1 - p.y_EN_ENT) * J1[N,rej] * sh.n_N_EN
    nothing
end

"""
Generalised reserve drain for any flux column *col* (ie gro)
and any combination of reserves.
"""
reserve_drain!(o, col, drain, θ) = begin
    J_CN = -drain * (1.0 - θ) # fraction on drain from C and N reserves
    o.J[C,col] = J_CN/o.shared.y_E_CH_NO
    o.J[N,col] = J_CN/o.shared.y_E_EN
    o.J[E,col] = -drain * θ
    nothing
end

"""
Generalised reserve loss to track carbon. 

Loss is distributed between general and C and N reserves by the fraction θE
"""
reserve_loss!(o, loss, θ) = begin
    p, v, sh = unpack(o)
    ee = loss * θ # fraction of loss from E reserve
    ecn = loss - ee # fraction on loss from C and N reserves
    ec = ecn/sh.y_E_CH_NO
    en = ecn/sh.y_E_EN
    o.J1[C,los] += ec + en + ee
    o.J1[N,los] += (ec, en, ee) ⋅ (sh.n_N_EC, sh.n_N_EN, sh.n_N_E)
    nothing
end

conversion_loss!(o, loss, θ, dest_ratio) = begin
    p, v, sh = unpack(o)
    ee = loss * θ # fraction of loss from E reserve
    ecn = loss - ee # fraction on loss from C and N reserves
    ec = ecn/sh.y_E_CH_NO
    en = ecn/sh.y_E_EN
    o.J1[C,los] += ec + en + ee - loss
    o.J1[N,los] += (ec, en, ee) ⋅ (sh.n_N_EC, sh.n_N_EN, (sh.n_N_E - 1/θ * dest_ratio))
end

# reserve_loss!(o, ::Type{Val{:C}}, loss) = o.J1[C,los] += loss
# reserve_loss!(o, ::Type{Val{:N}}, loss) = o.J1[N,los] += loss/o.shared

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

stoich_merge_losses(Ja, Jb, JEab, n_a, n_b, n_Eab) = begin 
    lossa = Ja + Jb - JEab
    lossb = Ja * n_a + Ja * n_b - JEab * n_Eab  
    lossa, lossb
end

"""
κtra is the difference paramsbetween κsoma and κmat
"""
κtra(o) = (1.0 - o.params.κsoma - κmat(o))

κmat(o::Organ) = κmat(o.params.maturity)
κmat(maturity::Maturity) = maturity.κmat
κmat(maturity::Nothing) = 0.0

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
feedback!(f::Nothing, o, u) = nothing
feedback!(f::Autophagy, o, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, rate(o.vars))
    aph = u[V] * (oneunit(hs) - hs)
    # TODO this should be a lossy process
    o.J[E,fbk] += aph
    o.J[V,fbk] -= aph
end

set_scaling!(o, u) = set_var!(o.vars, :scale, scaling(o.params.scaling, u[V]))

scaling(f::KooijmanArea, uV) = begin
    uV > zero(uV) || return 1.0 # || error("Mass is less than zero, I think its dead...")
    (uV / f.M_Vref)^(-uV / f.M_Vscaling)
end
scaling(f::Nothing, uV) = 1.0

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
germinated(M_V, M_Vgerm) = M_V > M_Vgerm


allometric_height(f::SqrtAllometry, o, u) = begin
    p, v, sh = unpack(o)
    dim = oneunit(u[V] * sh.w_V)
    sqrt((u[P] * sh.w_P + u[V] * sh.w_V) / dim) * f.allometry
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
# ┏━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃     ┃ catS       │ rejS      │ losS      ┃     ┃ catR       │ rejR      │ losR      ┃
# ┣━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃EES  ┃ J_EES,catS │ 0         │ 0         ┃EES  ┃ J_EER,catR │ 0         │ 0         ┃
# ┃CNS  ┃ J_CNS,catS │ 0         │ 0         ┃NS   ┃ J_CNR,catR │ 0         │ 0         ┃
# ┃CS   ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃CS   ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
# ┃NS   ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃ENS  ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
# ┃ES   ┃ J_ES,catS  │ 0         │ 0         ┃ES   ┃ J_ER,catR  │ 0         │ 0         ┃
# ┗━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
#
#
#
#
#
#
#

# M: State vector diagram.
# ┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
