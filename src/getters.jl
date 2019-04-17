# We use getters so the model never accesses struct fields directly, 
# except in methods for specific leaf types where declaring all the gettier 
# methods would be overkill.
#
# A benefit of this is parameters can be easily moved between shared and specific 
# parameter structs.


# So we don't have to depend on all of Lazy.jl
macro forward(ex, fs)
  @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
  T = esc(T)
  fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
  :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
       for f in fs]...);
    nothing)
end


@inline dead(o::AbstractOrganism) = o.dead[]
@inline set_dead!(o::AbstractOrganism, val) = o.dead[] = val
@inline environment(o::AbstractOrganism) = o.environment

@inline vars(o::AbstractOrgan) = o.vars
@inline flux(o::AbstractOrgan) = o.J
@inline flux1(o::AbstractOrgan) = o.J1

@inline tstep(v) = v.t[1]
@inline set_tstep!(v, val) = v.t[1] = val  

const varfuncs = [:shape, :rate, :temp, :θE, :tempcorrection, :height]
for (v, sv) in zip(varfuncs, Symbol.(:set_, varfuncs, :!))
    @eval @inline ($v)(vars) = vars.$v[tstep(vars)]
    @eval @inline ($sv)(vars, val) = vars.$v[tstep(vars)] = val  
end

@inline depth(v) = height(v)

rate_formula(p) = p.rate_formula
assimilation_pars(p) = p.assimilation_pars
shape_pars(p) = p.shape_pars
allometry_pars(p) = p.allometry_pars
maturity_pars(p) = p.maturity_pars
trans_pars(p) = p.trans_pars
production_pars(p) = p.production_pars
rejection_pars(p) = p.rejection_pars
germination_pars(p) = p.germination_pars
turnover_pars(p) = p.turnover_pars
feedback_pars(p) = p.feedback_pars
tempcorr_pars(p) = p.tempcorr_pars
core_pars(p) = p.core_pars
su_pars(p) = p.su_pars
catabolism_pars(p) = p.catabolism_pars
maintenance_pars(p) = p.maintenance_pars

@forward AbstractOrgan.vars θE, temp, set_temp!, tempcorrection, set_tempcorrection!, 
         height, set_height!, rate, set_rate!, shape, set_shape!, tstep, set_tstep!

@forward AbstractOrgan.params rate_formula, assimilation_pars, shape_pars, allometry_pars, maturity_pars,
                              trans_pars, production_pars, rejection_pars, germination_pars, turnover_pars

@forward AbstractOrgan.shared maintenance_pars, feedback_pars, su_pars, tempcorr_pars, catabolism_pars, core_pars,
                              y_V_E, y_E_EC, y_E_EN, n_N_P, n_N_V, n_N_E, n_N_EC, n_N_EN, w_V, w_C, w_N, w_E

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

j_E_mai(o::AbstractOrgan) = j_E_mai(maintenance_pars(o))

κtra(o::AbstractOrgan) = κtra(trans_pars(o))
κtra(o::Nothing) = 0.0

κmat(o::AbstractOrgan) = κmat(maturity_pars(o))
κmat(::Nothing) = 0.0

κsoma(o::AbstractOrgan) = oneunit(κtra(o)) - κtra(o) - κmat(o)

mass(o, u) = u.V * w_V(o)
