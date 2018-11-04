
Defaults.get_default(t::Type) = begin 
    d = default(t) 
    u = units(t)
    add_units.(d, u)
end

add_units(::Nothing, u) = nothing
add_units(x, ::Nothing) = x
add_units(::Nothing, ::Nothing) = nothing
add_units(x::Number, u::Unitful.Units) = x * u
add_units(x::AbstractArray, u::Unitful.Units) = x .* u

build_vars(vars, time) = begin
    len = length(time)
    len == length(vars.rate) && return vars 

    fields = []
    for fname in fieldnames(typeof(vars))
        ft = fieldtype(typeof(vars), fname)
        if ft <: AbstractArray
            push!(fields, fill(getfield(vars, fname)[1], len)) 
        else
            push!(fields, getfield(vars, fname))
        end
    end
    typeof(vars).name.wrapper(fields...)
end
build_J(one_flux, T, time) = build_flux(one_flux, T, STATE, TRANS, time)
build_J1(one_flux, T, time) = build_flux(one_flux, T, STATE1, TRANS1, time)
build_flux(one_flux, T, x, y, time) = zeros(typeof(one_flux), length(x), length(y), length(time))

split_state(o::Tuple, u::AbstractArray) = split_state(o, u, 0)
split_state(o::Tuple{O,Vararg}, u::AbstractArray, offset) where O = begin
    v = view(u, offset+1:offset+STATELEN)
    lv = LVector{eltype(v),typeof(v),STATE}(v)
    (lv, split_state(tail(o), u, offset + STATELEN)...)
end
split_state(o::Tuple{}, u::AbstractArray, offset) = ()

" Sum flux matrix " 
sum_flux!(du, o::Tuple) = begin
    offset_apply!(sum_flux!, du, o, 0)
    du
end
sum_flux!(du, organ::Organ, offset::Int) = begin
    J = organ.J
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

" update vars and flux to record for time t "
define_organs(o::Organism, t) = define_organs(o.params, o.records, o, t)
define_organs(params::Tuple{P,Vararg}, records::Tuple{R,Vararg}, organism, t) where {P,R} = begin
    rec = records[1]
    t = length(rec.vars.rate) != 1 ? floor(Int, ustrip(t)) + 1 : 1
    rec.vars.t = t
    vJ = view(rec.J, :, :, t)
    vJ1 = view(rec.J1, :, :, t)
    J = LMatrix{eltype(vJ),typeof(vJ),STATE,TRANS}(vJ)
    J1 = LMatrix{eltype(vJ1),typeof(vJ1),STATE1,TRANS1}(vJ1)
    organ = Organ(params[1], organism.shared, rec.vars, J, J1)
    (organ, define_organs(tail(params), tail(records), organism, t)...)
end
define_organs(params::Tuple{}, records::Tuple{}, organism, t) = ()

check_params(o::Organ) = begin 
    p = o.params; sh = o.shared
    p.y_P_V <= sh.n_N_P/sh.n_N_V || error("y_P_V too high for N conservation ", (p.y_P_V, sh.n_N_P, sh.n_N_V))
    p.y_V_E <= sh.n_N_V/sh.n_N_E || error("y_V_E too high for these valuse of n_N_V and n_N_E ", (p.y_V_E, sh.n_N_V, sh.n_N_E ))
    2 <= p.y_E_CH_NO + p.y_E_EN || error("y_ECH_NO or y_E_EN too high for C conservation ", (p.y_E_CH_NO, p.y_E_EN))
    2sh.n_N_E <= sh.n_N_EC * p.y_E_CH_NO + sh.n_N_EN * p.y_E_EN || error("y_ECH_NO or y_E_EN too high for N conservation ", 
                                                                         (2sh.n_N_E, sh.n_N_EC, p.y_E_CH_NO, sh.n_N_EN, p.y_E_EN))

    # These will be required if structures can have different reserve N ratios
    # p.y_E_ET < n_N_E/n_N_E
    # p.y_EC_ECT <n_N_C/n_N_C
    # p.y_EN_ENT < n_N_N/n_N_N
end

