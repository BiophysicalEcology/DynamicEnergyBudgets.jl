assimilation_vars(v) = v.assimilation_vars
shape(v) = v.shape[v.t]
rate(v) = v.rate[v.t]
temp(v) = v.temp[v.t]
θE(v::Vars) = v.θE[v.t]
θE(o::Organ) = θE(o.vars)
tempcorrection(v) = v.tempcorrection[v.t]
height(v) = v.height[v.t]
set_var!(o::Organ, fname, val) = set_var!(o.vars, fname, val)
set_var!(v, fname, val) = getfield(v, fname)[v.t] = val

name(o::Organ) = o.params.name
rate_formula(o::Organ) = o.params.rate_formula
assimilation_pars(o::Organ) = o.params.assimilation_pars
shape_pars(o::Organ) = o.params.shape_pars
allometry_pars(o::Organ) = o.params.allometry_pars
maturity_pars(o::Organ) = o.params.maturity_pars
trans_pars(o::Organ) = o.params.trans_pars
production_pars(o::Organ) = o.params.production_pars
rejection_pars(o::Organ) = o.params.rejection_pars
germination_pars(o::Organ) = o.params.germination_pars
turnover_pars(o) = o.params.turnover_pars

feedback_pars(o::Organ) = o.shared.feedback_pars
su_pars(o::Organ) = o.shared.su_pars
tempcorr_pars(o::Organ) = o.shared.tempcorr_pars
core_pars(o::Organ) = o.shared.core_pars

j_E_mai(o::Organ) = o.params.maintenance_pars.j_E_mai
j_P_mai(o::Organ) = production_pars(o).j_P_mai
y_P_V(o::Organ) = production_pars(o).y_P_V
M_Vgerm(o::Organ) = germination_pars(o).M_Vgerm
y_V_E(o::Organ) = core_pars(o).y_V_E
y_E_EC(o::Organ) = core_pars(o).y_E_EC
y_E_EN(o::Organ) = core_pars(o).y_E_EN
n_N_P(o::Organ) = production_pars(o).n_N_P
n_N_V(o::Organ) = core_pars(o).n_N_V
n_N_E(o::Organ) = core_pars(o).n_N_E
n_N_EC(o::Organ) = core_pars(o).n_N_EC
n_N_EN(o::Organ) = core_pars(o).n_N_EN
w_V(o::Organ) = core_pars(o).w_V
w_C(o::Organ) = core_pars(o).w_C
w_N(o::Organ) = core_pars(o).w_N
w_E(o::Organ) = core_pars(o).w_E

κtra(o::Organ) = κtra(trans_pars(o))
κtra(trans_pars::AbstractTranslocation) = trans_pars.κtra
κtra(o::Nothing) = 0.0

κmat(o::Organ) = κmat(maturity_pars(o))
κmat(maturity::Maturity) = maturity_pars.κmat
κmat(maturity::Nothing) = 0.0

κsoma(o::Organ) = (1.0 - κtra(o) - κmat(o))
