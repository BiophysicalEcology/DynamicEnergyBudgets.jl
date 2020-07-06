"""
    half_saturation(max, half, x)

Half satration curve.
"""
half_saturation(max, half, x) = max/(oneunit(half/x) + half/x)

"""
    reserve_drain!(o::AbstractOrgan, column, drain)

Generalised reserve drain for any flux column (ie :gro), 
specified by its `Symbol` or `Int` index.
"""
@inline reserve_drain!(o::AbstractOrgan, col, drain) = reserve_drain!(has_reserves(o), o, col, drain)
@inline reserve_drain!(::HasCNE, o, col, drain) = begin
    ΘE, J = θE(o), flux(o)
    J_CN = -drain * (1 - θ) # fraction on drain from C and N reserves
    @inbounds J[:C,col] = J_CN / y_E_C(o)
    @inbounds J[:N,col] = J_CN / y_E_N(o)
    @inbounds J[:E,col] = -drain * ΘE
    nothing
end
@inline reserve_drain!(::HasCN, o, col, drain) = begin
    J = flux(o)
    @inbounds J[:C,col] = -drain / y_E_C(o)
    @inbounds J[:N,col] = -drain / y_E_N(o)
    nothing
end

