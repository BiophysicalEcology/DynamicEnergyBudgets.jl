using Base: tail
import CompositeFieldVectors: reconstruct
using ForwardDiff: Dual

# Parameter helper functions
watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L 
fraction_per_litre_gas_to_mols(frac) = frac / 22.4

Defaults.get_default(t::Type) = begin 
    d = default(t) 
    u = units(t)
    combine.(u, d)
end

combine(::Void, ::Void) = nothing
combine(a, ::Void) = a
combine(::Void, b) = b
combine(a, b) = a * b


build_record(a, time) = [deepcopy(a) for i = 1:length(time)]
build_J(one_flux, T) = build_flux(one_flux, T, STATE, TRANS)
build_J1(one_flux, T) = build_flux(one_flux, T, STATE1, TRANS1)
build_flux(one_flux, T, x, y) = T[zero(one_flux) for a in 1:length(x), b in 1:length(y)]

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
    i = length(rec.vars) > 1 ? floor(Int, ustrip(t)) + 1 : 1
    organ = Organ(params[1], organism.shared, rec.vars[i], rec.J[i], rec.J1[i])
    (organ, define_organs(tail(params), tail(records), organism, t)...)
end
define_organs(params::Tuple{}, records::Tuple{}, organism, t) = ()


(o::Organism)(du, u, p::Void, t::Number) = o(du, u, t, define_organs(o, t))
(o::Organism)(du::AbstractVector{<:Dual}, u::AbstractVector{<:Dual}, p::AbstractVector{<:Real}, t::Number) = begin
    # println("p not Duuaaaaalllllll")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = du .* u"mol/hr"
    rec = dualize_records(o, du1, u)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    du .= ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Dual}, u::AbstractVector{<:Dual}, p::AbstractVector{<:Dual}, t::Number) = begin
    # println("Alllllll Duuaaaaallllll")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = [zero(p[1]) .* u"mol/hr" for d in du]
    rec = dualize_records(o, du1, u)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    du .= ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Dual}, t::Number) = begin
    # println("Reset duuuuuuuuuu")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = [zero(p[1]) .* u"mol/hr" for d in du]
    rec = dualize_records(o, du1, p)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Number) = begin
    # println("Parammmmmmmmmmmmmmmmms")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = du .* u"mol/hr"
    rec = Records(o.params)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    du2 = ustrip.(du1)
    if eltype(du2) == eltype(du)
        du .= du2
    end
    du2
end

dualize_records(o, du, u) = begin
    vars1 = dualize_vars(o.records[1].vars[1], zero(u[1]))
    vars2 = dualize_vars(o.records[2].vars[1], zero(u[1]))
    Records(o.params, vars=(vars1, vars2), val=zero(du[1]))
end

dualize_vars(vars, val) = begin
    flatvars = flatten(vars)
    dualvars = ForwardDiff.seed!(fill(val, length(flatvars)), flatvars)
    reconstruct(vars, dualvars)
end

