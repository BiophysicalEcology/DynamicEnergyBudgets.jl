"""
Run a DEB organism model.
"""
(o::Organism)(du, u, p::Void, t::Number) = o(du, u, t, define_organs(o, t))

function (o::Organism)(du, u, t::Number, organs)
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
function debmodel!(organs::Tuple, u::Tuple)
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
    dissipation!(o, u)
    feedback!(o.shared.feedback, o, u)
end

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not finalise flux in J - operates only on J1 (intermediate storage)
"""
function catabolism!(o, u)
    p, v, sh, J, J1 = unpack(o); va = assimilation(v);

    turnover = (p.k_EC, p.k_EN, p.k_E) .* tempcorrection(v) .* scale(v)
    m = (u[C], u[N], u[E]) ./ u[V]
    j_E_mai = p.j_E_mai * tempcorrection(v)

    rate = find_rate(m, turnover, j_E_mai, sh.y_E_CH_NO, sh.y_E_EN, p.y_V_E, p.κsoma)
    set_var!(v, :rate, rate)
    # rate > zero(rate) || error("rate is less than zero, thermodyamics disagrees")

    (J1[C,cat], J1[N,cat], J1[EE,cat]) = catabolic_flux((u[C], u[N], u[E]), turnover, rate)
    (J1[C,rej], J1[N,rej], J1[CN,cat], o.J1[C,los], o.J1[N,los]) =
        stoich_merge(J1[C,cat], J1[N,cat], sh.y_E_CH_NO, sh.y_E_EN)

    J1[E,cat] = J1[EE,cat] + J1[CN,cat] # Total catabolic flux
    θE(v) = J1[EE,cat]/J1[E,cat] # Proportion of general reserve flux in total catabolic flux
    nothing
end

"""
Dissipation for any reserve.
Growth, maturity and maintenence are grouped as dissipative processes.
"""
function dissipation!(o, u)
    growth!(o, u)
    maturity!(o, u)
    maintenence!(o, u)
    product!(o, u)
end

"""
Allocates reserves to growth.
"""
function growth!(o, u)
    p, v, sh, J, J1 = unpack(o)
    grow = rate(v) * u[V]
    J[V,gro] = grow 
    drain = -(1/p.y_V_E) * grow 
    loss = (1/p.y_V_E - 1) * rate(v) * u[V]
    reserve_drain!(o, gro, drain, θE(v))
    reserve_loss!(o, loss)
end

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(o.params.maturity, o, u)
maturity!(f::Maturity, o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    J[M,gro] = f.κmat * J1[E,cat]
    maint = -f.j_E_mat_mai * tempcorrection(v) * u[V]
    drain = -J[M,gro] + maint # min(u[V], f.M_Vmat))
    reserve_drain!(o, mat, drain, θE(v))
    reserve_loss!(o, -maint)
end
maturity!(f::Void, o, u) = nothing

"""
Allocates reserve drain due to maintenance.
"""
function maintenence!(o, u)
    drain = -o.params.j_E_mai * tempcorrection(o.vars) * u[V]
    reserve_drain!(o, mai, drain, θE(o.vars))
    reserve_loss!(o, -drain) # all maintenance is loss
end

"""
Allocates waste products from growth and maintenance.
"""
function product!(o, u)
    p, v, sh, J, J1 = unpack(o)
    J[P,gro] = J[V,gro] * p.y_P_V
    J[P,mai] = u[V] * p.j_P_mai * tempcorrection(v)
    # undo the reserve loss from growth: it went to product
    loss_correction = (J[P,gro] + J[P,mai])
    reserve_loss!(o, loss_correction)
end


"""
Translocation occurs between adjacent organs. 
This function is identical both directiono, so on represents
whichever is not the current organs. 

Will not run with less than 2 organs.
"""
translocation!(organs::Tuple{}) = nothing
translocation!(organs::Tuple{Organ}) = nothing
translocation!(organs::Tuple{Organ, Organ}) = begin
    reuse_rejected!(organs[1], organs[2], 1.0)
    reuse_rejected!(organs[2], organs[1], 1.0)
    translocate!(organs[1], organs[2], 1.0)
    translocate!(organs[2], organs[1], 1.0)
end
# translocation!(organs::Tuple) = translocation!(organs, organs)

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
buildprops(o::Organ) = buildprops(o.params.translocation.proportions)
buildprops(x::Void) = (1.0)
buildprops(x::Number) = (x, 1 - x)
buildprops(xs::Tuple) = (xs..., 1 - sum(xs))

"""
Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs. 
This function is identical both directiono, and ox represents
whichever is not the current organs. Will not run with less than 2 organs.
"""
function translocate!(source, dest, prop)
    trans = κtra(source) * prop * source.J1[E,cat]
    loss = (1 - source.params.y_E_ET) * trans
    reserve_drain!(source, tra, -trans, θE(source.vars))
    reserve_loss!(source, loss)
    # incoming translocation
    transx = κtra(dest) * source.J1[E,cat]
    source.J[E,tra] += dest.params.y_E_ET * transx
    nothing
end

"""
Reallocate state rejected from synthesizing units.

TODO add a 1-organs method
Also how does this interact with assimilation?
"""
function reuse_rejected!(source, dest, prop)
    p = source.params
    # rejected reserves are translocated and used in assimilation.
    if typeof(source.params.assimilation) <: AbstractCAssim 
        # source.J[C,rej] = -source.J1[C,rej]
        # dest.J[C,tra] = p.y_EC_ECT * source.J1[C,rej]
        # source.J[N,rej] = source.J1[C,rej]
    elseif typeof(source.params.assimilation) <: AbstractNAssim 
        source.J[N,rej] = -source.J1[N,rej]
        dest.J[N,tra] = p.y_EN_ENT * source.J1[N,rej]
        # source.J[C,rej] = source.J1[C,rej]
    end
    source.J1[C,los] += (1 - p.y_EC_ECT) * source.J1[C,rej]
    source.J1[N,los] += (1 - p.y_EN_ENT) * source.J1[N,rej]
    nothing
end

"""
Generalised reserve drain for any flux column *col* (ie gro)
and any combination of reserves.
"""
function reserve_drain!(o, col, drain, θE)
    J_CN = drain * (1.0 - θE) # fraction on drain from C and N reserves
    o.J[C,col] = J_CN/o.shared.y_E_CH_NO
    o.J[N,col] = J_CN/o.shared.y_E_EN
    o.J[E,col] = drain * θE
    nothing
end

"""
Generalised reserve loss to track carbon. 
"""
function reserve_loss!(o, loss)
    o.J1[E,los] += loss
    nothing
end

reserve_loss!(o, ::Type{Val{:C}}, loss) = o.J1[C,los] += loss
# reserve_loss!(o, ::Type{Val{:N}}, loss) = o.J1[N,los] += loss/o.shared

"""
    stoich_merge(Ja, Jb, ya, yb)
Merge fluxes stoichiometrically into general reserve Eab based on yeild 
fractions ya and yb. An unmixed proportion is returned as unmixed reserves Ea and Eb.
Losses are also calculated.
"""
stoich_merge(Ja, Jb, ya, yb) = begin
    JEab = synthesizing_unit(Ja * ya, Jb * yb) 
    JEa = Ja - JEab/ya                        
    JEb = Jb - JEab/yb                     
    lossa = Ja + Jb - JEa - JEb - JEab
    lossb = zero(lossa)
    (JEa, JEb, JEab, lossa, lossb)
end

"""
κtra is the difference paramsbetween κsoma and κmat
"""

κtra(o) = (1.0 - o.params.κsoma - κmat(o)) 

κmat(o::Organ) = κmat(o.params.maturity)
κmat(maturity::Maturity) = maturity.κmat
κmat(maturity::Void) = 0.0

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

Without a function like this you will likely be occasionally breaking the 
laws of thermodynamics by introducing negative rates.
"""
feedback!(f::Void, o, u) = nothing
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
scaling(f::Void, uV) = 1.0

"""
Check if germination has happened. Independent for each organ,
although this may not make sense. A curve could be better for this too.
"""
germinated(M_V, M_Vgerm) = M_V > M_Vgerm 


"""
"""
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
# ┏━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃     ┃ catS       │ rejS      │ losS      ┃
# ┣━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃EES  ┃ J_EES,catS │ 0         │ 0         ┃
# ┃CNS  ┃ J_CNS,catS │ 0         │ 0         ┃
# ┃CS   ┃ J_CS,catS  │ J_CS,rejS │ J_CS,losS ┃
# ┃NS   ┃ J_NS,catS  │ J_NS,rejS │ J_NS,losS ┃
# ┃ES   ┃ J_ES,catS  │ 0         │ 0         ┃
# ┃EES  ┃ J_EER,catR │ 0         │ 0         ┃
# ┃NS   ┃ J_CNR,catR │ 0         │ 0         ┃
# ┃CS   ┃ J_ECR,catR │ J_CR,rejR │ J_CR,losR ┃
# ┃ENS  ┃ J_ENR,catR │ J_NR,rejR │ J_NR,losR ┃
# ┃ES   ┃ J_ER,catR  │ 0         │ 0         ┃
# ┗━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# M: State vector diagram.
# ┏━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃                      ┃ MS SubArray            ┃ MR SubArray            ┃
# ┃State variable (mols) ┃ PS │ VS │ CS │ NS │ ES ┃ PR │ VR │ CR │ NR │ ER ┃
# ┗━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┛
