import Flatten
using ForwardDiff: Dual

dualize_records(o, du, u) = begin
    vars1 = dualize_vars(o.records[1].vars[1], zero(u[1]))
    vars2 = dualize_vars(o.records[2].vars[1], zero(u[1]))
    Records(o.params, vars=(vars1, vars2), val=zero(du[1]))
end

dualize_vars(vars, val) = begin
    flatvars = flatten(Vector, vars)
    dualvars = ForwardDiff.seed!(fill(val, length(flatvars)), flatvars)
    reconstruct(vars, dualvars)
end

(o::Organism)(du::AbstractVector{<:Dual}, u::AbstractVector{<:Dual}, p::AbstractVector{<:Real}, t::Number) = begin
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = du .* u"mol/hr"
    rec = dualize_records(o, du1, u)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    du .= ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Dual}, u::AbstractVector{<:Dual}, p::AbstractVector{<:Dual}, t::Number) = begin
    println("Dual/Dual")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = [zero(p[1]) .* u"mol/hr" for d in du]
    rec = dualize_records(o, du1, u)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    du .= ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Dual}, t::Number) = begin
    println("Real/Dual")
    u1 = u .* u"mol"
    o1 = reconstruct(o, p)
    du1 = [zero(p[1]) .* u"mol/hr" for d in du]
    rec = dualize_records(o, du1, p)
    o1(du1, u1, t * u"hr", define_organs(o1.params, rec, o1, t))
    ustrip.(du1)
end
(o::Organism)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Number) = begin
    println("Real/Real")
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

function deb_function(sol)
   tot_loss = 0.0
   if any((s.retcode != :Success for s in sol))
     tot_loss = Inf
   else
     # calculation for the loss here
   end
   tot_loss
end

