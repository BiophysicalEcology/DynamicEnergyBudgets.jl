import CompositeFieldVectors: reconstruct

"""
Run a DEB organism model.
"""
(o::Organism)(du, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Number; record=false) = begin
    Handle unitless inputs
    du1 = Any[d .* u"mol/hr" for d in du]
    u1 = u .* u"mol"
    reconstruct(o, p)(du1, u1, t * u"hr") 
    if typeof(du1[1]) <: ForwardDiff.Dual
        ustrip.(du1)
    else
        du2 = ustrip.(du1)
        if eltype(du2) == eltype(du)
            du .= du2
        end
        du2
    end
end
(o::Organism)(du, u, p::Void, t::Number; record=false;) = begin
    record && keep_records!(o, t)
    o(du, u, t) 
end

function (o::Organism)(du, u, t::Number)
    ux = split_state(o, u)
    organs = o.organs
    # apply_environment!(o, ux, t)
    apply_environment!(organs[1].params.assimilation, organs[1], u, o.environment, t)
    apply_environment!(organs[2].params.assimilation, organs[2], u, o.environment, t)
    debmodel!(o.organs, ux)
    sum_flux!(du, o)
end

"""
A generalised multi-reserve, multi-organ Dynamic Energy Budget model.

Applies metabolism, translocation and assimilation mehtods to N organs.

settings is a struct with required model data, DEBSettings or similar.
t is the timestep
"""
function debmodel!(organs::Tuple, u::Tuple)
    swapped = (Base.tail(organs)..., organs[1])
    metabolism!(organs, u)
    run_translocation!(organs, swapped)
    assimilation!(organs, swapped, u)
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
    feedback!(o, o.shared.feedback, u)
end

"""
    catabolism!(o, u, t::Number)
Catabolism for E, C and N, or C, N and E reserves.
Does not finalise flux in J - operates only on J1 (intermediate storage)
"""
function catabolism!(o, u)
    p, v, sh, J, J1 = unpack(o); va = v.assimilation;
    rate = v.rate
    tempscaling = v.tempcorr .* v.scale

    scaled_turnover = (p.k_EC, p.k_EN, p.k_E) .* tempscaling
    m = (u[C], u[N], u[E]) ./ u[V]
    j_E_mai = p.j_E_mai * v.tempcorr
    v.rate = find_rate(m, scaled_turnover, j_E_mai, sh.y_E_CH_NO, sh.y_E_EN, p.y_V_E, p.κsoma)

    (J1[C,cat], J1[N,cat], J1[EE,cat]) = 
        catabolic_fluxes((u[C], u[N], u[E]), scaled_turnover, rate)
    (J1[C,rej], J1[N,rej], J1[CN,cat]) =
        synthesizing_unit(J1[C,cat], J1[N,cat], sh.y_E_CH_NO, sh.y_E_EN)

    J1[E,cat] = J1[EE,cat] + J1[CN,cat] # Total catabolic flux
    v.θE = J1[EE,cat]/J1[E,cat] # Proportion of general reserve flux in total catabolic flux
    nothing
end

"""
Dissipation for any reserve.
Growth, maturity and maintenence are grouped as dissipative processes.
"""
function dissipation!(o, u)
    growth!(o, u)
    maturity!(o.params.maturity, o, u)
    maintenence!(o, u)
    product!(o, u)
end

"""
Allocates reserves to growth.
"""
function growth!(o, u)
    p, v, sh, J, J1 = unpack(o)
    grow = o.vars.rate * u[V]
    J[V,gro] = grow 
    drain = -(1/p.y_V_E) * grow 
    loss = (1/p.y_V_E - 1) * v.rate * u[V]
    reserve_drain!(o, gro, drain, v.θE)
    reserve_loss!(o, loss)
end

"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(f::Maturity, o, u) = begin
    p, v, sh, J, J1 = unpack(o)
    J[M,gro] = f.κrep * J1[E,cat]
    maint = -f.j_E_rep_mai * v.tempcorr * u[V]
    drain = -J[M,gro] + maint # min(u[V], f.M_Vrep))
    reserve_drain!(o, rep, drain, v.θE)
    reserve_loss!(o, -maint)
end

# function maturity!(f, o, u) end
# maturity!(f::Maturity, o, u::U) where U = maturity!(f, o, u, has_M(U))
# maturity!(f::Maturity, o, u, ::NoM) = begin
#     p, v, sh, J, J1 = unpack(o)
#     # TODO: why does rep maintenance stop increasing at M_Vrep?
#     # Is this a half finished reproduction model?
#     drain = -(f.κrep * J1[E,cat] + -f.j_E_rep_mai * v.tempcorr * min(u[V], f.M_Vrep))
#     reserve_drain!(o, rep, drain, v.θE)
#     reserve_loss!(o, -drain)
# end

"""
Allocates reserve drain due to maintenance.
"""
function maintenence!(o, u)
    drain = -o.params.j_E_mai * o.vars.tempcorr * u[V]
    reserve_drain!(o, mai, drain, o.vars.θE)
    reserve_loss!(o, -drain) # all maintenance is loss
end

"""
Allocates waste products from growth and maintenance.
"""
function product!(o, u)
    p, v, sh, J, J1 = unpack(o)
    J[P,gro] = J[V,gro] * p.y_P_V
    J[P,mai] = u[V] * p.j_P_mai * v.tempcorr
    # undo the reserve loss from growth: it went to product
    loss_correction = -(J[P,gro] + J[P,mai])
    reserve_loss!(o, loss_correction)
end


"""
Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs. 
This function is identical both directiono, so on represents
whichever is not the current organs. Will not run with less than 2 organs.

FIXME this will be broken for organs > 2
"""

translocation!(organs::Tuple{}, swapped) = nothing
translocation!(organs::NTuple{1}, swapped) = nothing
translocation!(organs::Tuple, swapped) = apply(translocation!, organs, swapped)
translocation!(o::Organ, on::Organ) = begin
    reuse_rejected!(o, on)
    translocate!(o, on)
end


"""
Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs. 
This function is identical both directiono, so on represents
whichever is not the current organs. Will not run with less than 2 organs.

FIXME this will be broken for organs > 2
"""
function translocate!(o, on)
    p = o.params
    trans = κtra(p) * o.J1[E,cat]
    loss = (1 - p.y_E_ET) * trans
    reserve_drain!(o, tra, -trans, o.vars.θE)
    reserve_loss!(o, loss)
    # incoming translocation
    transn = κtra(on.params) * o.J1[E,cat]
    o.J[E,tra] += on.params.y_E_ET * transn
    nothing
end

"""
Reallocate state rejected from synthesizing units.

TODO add a 1-organs method
Also how does this interact with assimilation?
"""
function reuse_rejected!(o, on)
    p = o.params
    # rejected reserves are translocated and used in assimilation.
    o.J[C,rej] = -o.J1[C,rej]
    o.J[N,rej] = -o.J1[N,rej]
    on.J[C,tra] = p.y_EC_ECT * o.J1[C,rej]
    on.J[N,tra] = p.y_EN_ENT * o.J1[N,rej]
    o.J1[C,los] += (1 - p.y_EC_ECT) * o.J1[C,rej]
    o.J1[N,los] += (1 - p.y_EN_ENT) * o.J1[N,rej]
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
    o.J1[C,los] += loss/o.shared.y_E_CH_NO
    o.J1[N,los] += loss/o.shared.y_E_EN
    nothing
end

"""
κtra is the difference paramsbetween κsoma and κrep
"""
κtra(params) = 1.0 - params.κsoma - κrep(params)
κrep(params) = params.maturity.κrep

"""
Calculate rate formula. TODO: use Roots.jl for this
"""
function find_rate(m, scaled_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)
    # bounds = (-10.0oneunit(v.rate), 20.0oneunit(v.rate)) #rate_window(args...) 

    # Find the type of rate so that diff and units work with the guess and atol
    one_r = oneunit(eltype(scaled_turnover))
    # Get rate with a zero finder
    let m=m, scaled_turnover=scaled_turnover, j_E_mai=j_E_mai, y_E_CH_NO=y_E_CH_NO, y_E_EN=y_E_EN, y_V_E=y_V_E, κsoma=κsoma
        find_zero((-2one_r, 1one_r), Secant(), one_r*1e-10, 100) do x
            rate_formula(x, m, scaled_turnover, j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)
        end
    end
end

"""
Function to apply feedback on growth the process, such as autopagy in resource shortage.

Without a function like this you will likely be occasionally breaking the 
laws of thermodynamics by introducing negative rates.
"""
feedback!(o, f::Void, u) = nothing
feedback!(o, f::Autophagy, u) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, o.vars.rate)
    autophagy = u[V] * (oneunit(hs) - hs)
    o.J[C,gro] += autophagy/o.shared.y_E_CH_NO
    o.J[N,gro] += autophagy/o.shared.y_E_EN
    o.J[V,gro] -= autophagy
    nothing
end


set_scaling!(o, u) = begin
    o.vars.scale = scaling(o.params.scaling, u[V])
end

"""
"""
scaling(f::KooijmanArea, uV) = begin
    uV > zero(uV) || return 0.0
    (uV / f.M_Vref)^(-uV / f.M_Vscaling)
end

scaling(f::Void, uV) = 1.0

"""
Check if germination has happened. Independent for each organ,
although this may not make sense.
"""
germinated(M_V, M_Vgerm) = M_V > M_Vgerm 


"""
"""
allometric_height(f::SqrtAllometry, o, u) = begin
    p, v, sh = unpack(o)
    dim = oneunit(u[V] * sh.w_V)
    sqrt((u[P] * sh.w_P + u[V] * sh.w_V) / dim) * f.size
end
    
unpack(o::Organ) = o.params, o.vars, o.shared, o.J, o.J1

# J: Flux matrix diagram.
# Rows: state.
# Columns: transformations
# ┏━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃    ┃ assS       │ groS       │ maiS       │ repS       │ rejS       │ traS       ┃ assR       │ groR       │ maiR       │ repR │ rejR       │ traR       ┃
# ┣━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
# ┃    ┃ JSS SubArray                                                                ┃ JRR SubArray                                                          ┃
# ┃    ┃                                                                             ┃                                                                       ┃
# ┃PS  ┃ 0          │ J_PS,groS  │ J_PS,maiS  │ 0          │ 0          │ 0          ┃ 0          │ J_PR,groR  │ J_PR,maiR  │ 0    │ 0          │ 0          ┃
# ┃VS  ┃ 0          │ J_VS,groS  │ 0          │ 0          │ 0          │ 0          ┃ 0          │ J_VR,groR  │ 0          │ 0    │ 0          │ 0          ┃
# ┃RS  ┃ 0          │ 0          │ 0          │ J_MS,groS  │ 0          │ 0          ┃ 0          │ 0          │ 0          │ 0    │ 0          │ 0          ┃
# ┃ECS ┃ J_ECS,assS │ J_ECS,groS │ J_ECS,maiS │ J_ECS,repS │ J_ECS,rejS │ J_ECS,traS ┃ J_ECR,assR │ J_ECR,groR │ J_ECR,maiR │ 0    │ J_ECS,rejR │ J_ECR,traR ┃
# ┃ENS ┃ J_ENS,assS │ J_ENS,groS │ J_ENS,maiS │ J_ENS,repS │ J_ENS,rejS │ J_ENS,traS ┃ J_ENR,assR │ J_ENR,groR │ J_ENR,maiR │ 0    │ J_ENS,rejR │ J_ENR,traR ┃
# ┃ES  ┃ J_ES,assS  │ J_ES,groS  │ J_ES,maiS  │ J_ES,repS  │ 0          │ J_ES,traS  ┃ J_ER,assR  │ J_ER,groR  │ J_ER,maiR  │ 0    │ 0          │ J_ER,traR  ┃
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
