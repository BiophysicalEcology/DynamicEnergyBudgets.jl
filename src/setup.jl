
FieldDefaults.get_default(t::Type) = begin 
    d = default(t) 
    u = units(t)
    add_units.(d, u)
end

add_units(::Nothing, u) = nothing
add_units(x, ::Nothing) = x
add_units(::Nothing, ::Nothing) = nothing
add_units(x::Number, u::Unitful.Units) = x * u
add_units(x::AbstractArray, u::Unitful.Units) = x .* u


split_state(o::Tuple, u::AbstractArray) = split_state(o, u, 0)
split_state(o::Tuple{O,Vararg}, u::AbstractArray, offset) where O = begin
    v = view(u, offset+1:offset+STATELEN)
    lv = LArray{STATE}(v)
    (lv, split_state(tail(o), u, offset + STATELEN)...)
end
split_state(o::Tuple{}, u::AbstractArray, offset) = ()


" Sum flux matrix " 
sum_flux!(du, organs::Tuple) = begin
    offset_apply!(sum_flux!, du, organs, 0)
    du
end
sum_flux!(du, organ::Organ, offset::Int) = begin
    J = flux(organ)
    z = zero(J[1,1])
    for i in 1:size(J, 1) 
        s = z 
        for j in 1:size(J, 2)
            s += J[i,j]
        end
        du[i+offset] = s
    end
    offset + size(J, 1)
end


update_tstep!(o::Organ, t) = set_tstep!(o, calc_tstep(o, t))

calc_tstep(o::Organ, t) = calc_tstep(vars(o), t)
calc_tstep(vars, t) = length(vars.rate) != 1 ? floor(Int, ustrip(t)) : 1

calc_envtime(o, t) = t + o.environment_start[]

zero_flux!(o) = begin
    o.J .= zero(eltype(o.J))
    o.J1 .= zero(eltype(o.J1))
end


(o::Plant)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Number) = begin
    println("Real/Real")
    u1 = u .* mol
    o1 = reconstruct(o, p)
    du1 = du .* (mol/hr)
    rec = Records(o.params)
    o1(du1, u1, t * (hr), define_organs(o1.params, rec, o1, t))
    du2 = ustrip.(du1)
    if eltype(du2) == eltype(du)
        du .= du2
    end
    du2
end
(o::Plant)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::Nothing, t::Number) = begin 
    t = t * (hr)
    du1 = du .* (mol/hr)
    u1 = u .* (mol)
    o(du1, u1, t, define_organs(o, t))
    du .= ustrip.(du1)
end
(o::Plant)(du, u, p::Nothing, t::Number) = o(du, u, t, define_organs(o, t))


check_params(o::Tuple) = apply(check_params, o)
check_params(o::Organ) = begin 
    y_V_E(o) <= n_N_V(o)/n_N_E(o) || error("y_V_E too high for these valuse of n_N_V and n_N_E ", (y_V_E(o), n_N_V(o), n_N_E(o) ))
    2 <= y_E_EC(o) + y_E_EN(o) || error("y_ECH_NO or y_E_EN too high for C conservation ", (y_E_EC(o), y_E_EN(o)))
    # 2n_N_E(o) <= n_N_EC(o) * y_E_EC(o) + n_N_EN(o) * y_E_EN(o) || error("y_ECH_NO or y_E_EN too high for N conservation ", 
    # (2n_N_E(o), n_N_EC(o), y_E_EC(o), n_N_EN(o), y_E_EN(o)))
    2n_N_E(o) <= y_E_EC(o) + y_E_EN(o) || error("y_ECH_NO or y_E_EN too high for N conservation ", 
                                                   (n_N_E(o), y_E_EC(o), y_E_EN(o)))
    check_params(production_pars(o), o::Organ)

    # These will be required if structures can have different reserve N ratios
    # y_E_ET(o) < n_N_E/n_N_E
    # y_EC_ECT(o) <n_N_C/n_N_C
    # y_EN_ENT(o) < n_N_N/n_N_N
end
check_params(p::Nothing, o::Organ) = nothing 
check_params(p::Production, o::Organ) = begin 
    p.y_P_V <= n_N_P(o)/n_N_V(o) || error("y_P_V too high for N conservation ", (p.y_P_V, p.n_N_P, n_N_V(o)))
    p.j_P_mai <= j_E_mai(o) || error("j_P_mai must be lower than j_E_mai ", (p.j_P_mai, j_E_mai(o)))
end
