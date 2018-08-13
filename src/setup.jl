using Base: tail

# Parameter helper functions
watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L 
fraction_per_litre_gas_to_mols(frac) = frac / 22.4

Defaults.get_default(t::Type) = begin 
    d = default(t) 
    u = units(t)
    add_units.(d, u)
end

add_units(::Void, u) = nothing
add_units(x, ::Void) = x
add_units(::Void, ::Void) = nothing
add_units(x::Number, u::Unitful.Units) = x * u
add_units(x::AbstractArray, u::Unitful.Units) = x .* u

build_vars(vars, time) = begin
    len = length(time)
    len == length(vars.rate) && return vars 

    fields = []
    for fname in fieldnames(vars)
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
build_flux(one_flux, y, time) = begin
    a = fill(one_flux, length(x), length(y), time)
    AxisArray(a, Axis{:state}(x), Axis{:transformations}(y), Axis{:time}(time))
end

split_state(o::Tuple, u::AbstractArray) = split_state(o, u, 0)
split_state(o::Tuple{O,Vararg}, u::AbstractArray, offset) where O = begin
    v = view(u, offset+1:offset+STATELEN)
    (v, split_state(tail(o), u, offset + STATELEN)...)
end
split_state(o::Tuple{}, u::AbstractArray, offset) = ()

" sum flux matrix " 
sum_flux!(du, o) = begin
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
    organ = Organ(params[1], organism.shared, rec.vars, view(rec.J, :, :,t), view(rec.J1, :, :, t))
    (organ, define_organs(tail(params), tail(records), organism, t)...)
end
define_organs(params::Tuple{}, records::Tuple{}, organism, t) = ()
