
"""
    split_state(o::Tuple, u::AbstractArray)

Split state vector into multiple views for each organ.

We do this recursively for type stability.
"""
split_state(o::Tuple, u::AbstractArray) = split_state(o, u, 0)
split_state(o::Tuple{O,Vararg}, u::AbstractArray, offset) where O = begin
    statedim = dims(flux(o[1]), X)
    v = view(parent(u), offset+1:offset+length(statedim))
    lv = DimensionalArray(v, (statedim,))
    (lv, split_state(tail(o), u, offset + length(statedim))...)
end
split_state(o::Tuple{}, u::AbstractArray, offset) = ()

unwrap(::Val{X}) where X = X
unwrap(::Type{Val{X}}) where X = X

update_tstep!(o::Organ, t) = set_tstep!(o, calc_tstep(o, t))

calc_tstep(o::Organ, t) = calc_tstep(vars(o), t)
calc_tstep(vars, t) = length(rate(vars)) != 1 ? floor(Int, ustrip(t)) : 1

calc_envtime(o, t) = t + o.environment_start[]

zero_flux!(o) = flux(o) .= zero(eltype(flux(o)))

"""
    check_params(o)

Check DEB parameters to avoid breaking mass balance. Can be pased an organ or
tuple of organs.
"""
function check_params end

check_params(organs::Tuple) = map(check_params, organs)
check_params(o::Organ) = begin 
    y_V_E(o) <= n_N_V(o)/n_N_E(o) || error("y_V_E too high for these valuse of n_N_V and n_N_E ", (y_V_E(o), n_N_V(o), n_N_E(o) ))
    1 >= y_E_C(o) || error("y_E_C must be less than or equal to 1: ", y_E_C(o))
    # 1 <= y_E_N(o) || error("y_E_N must be less than or equal to 1: ", y_E_N(o))
    # 2n_N_E(o) <= n_N_C(o) * y_E_C(o) + n_N_N(o) * y_E_N(o) || error("y_E_C or y_E_N too high for N conservation ", 
    # (2n_N_E(o), n_N_C(o), y_E_C(o), n_N_N(o), y_E_N(o)))
    2n_N_E(o) <= y_E_C(o) + y_E_N(o) || error("y_E_C or y_E_N too high for N conservation ", 
                                              (n_N_E(o), y_E_C(o), y_E_N(o)))
    check_params(production_pars(o), o)

    # These will be required if structures can have different reserve N ratios
    # y_E_ET(o) < n_N_E/n_N_E
    # y_E_ECT(o) <n_N_C/n_N_C
    # y_N_NT(o) < n_N_N/n_N_N
end
check_params(p::Nothing, o::Organ) = nothing 
check_params(p::Production, o::Organ) = begin 
    p.y_P_V <= p.n_N_P/n_N_V(o) || error("y_P_V too high for N conservation ", (p.y_P_V, p.n_N_P, n_N_V(o)))
    p.j_P_mai <= j_E_mai(o) || error("j_P_mai must be lower than j_E_mai ", (p.j_P_mai, j_E_mai(o)))
end


""" 
    sum_flux!(du, organs::Tuple)

Sum flux matrix and write to `du`.
""" 
sum_flux!(du, organs::Tuple) = begin
    offset_apply(sum_flux!, du, organs, 0)
    du
end
sum_flux!(du, organ::Organ, offset::Int) = begin
    J = flux(organ)
    z = zero(J[1, 1])
    for i in 1:size(J, 1) 
        s = z 
        for j in 1:size(J, 2)
            s += J[i, j]
        end
        du[i + offset] = s
    end
    offset + size(J, 1)
end

using Base: tail

offset_apply(f, a::AbstractArray, o::Tuple{O,Vararg}, offset::Int, args...) where O = begin
    offset = f(a, o[1], offset, args...)
    offset_apply(f, a, tail(o), offset::Int, args...)
    nothing
end
offset_apply(f, a::AbstractArray, o::Tuple{}, offset::Int, args...) = nothing
