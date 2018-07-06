oneunit_flux(params, state) = oneunit(state[1] * params.k_E)

build_J(one_flux, state) = build_flux(one_flux, fieldnames(state), TRANS)
build_J1(one_flux, state) = build_flux(one_flux, fieldnames(state), TRANS1)
build_flux(one_flux, x, y) = Array(fill(zero(one_flux), length(x), length(y)))
# build_flux(one_flux, x, y) = begin
#     a = Any[0.0 * one_flux for a in 1:length(x), b in 1:length(y)]
#     AxisArray(a, Axis{:state}(x), Axis{:transformations}(y))
# end

build_record(a, time) = 
    AxisArray([deepcopy(a) for i = 1:length(time)], Axis{:time}(time))

" copy the diffeq state to organs "

set_state!(o::Organism, u) = offset_apply!(set_state!, o.organs, u, 0)
set_state!(organ::Organ, u::AbstractArray, offset::Int) = begin
    typ = typeof(organ.state).name.wrapper
    organ.state = typ(u[1+offset:length(organ.state)+offset])
    offset + length(organ.state)
end

split_state(o::Organism, u::AbstractArray) = split_state(o.organs, u, 0)
split_state(o::Tuple{O,Vararg}, u::AbstractArray, offset) where O = 
    (view(u, offset+1:offset+STATELEN), split_state(Base.tail(o), u, offset + STATELEN)...)
split_state(o::Tuple{}, u::AbstractArray, offset) = tuple()

" sum flux matrix " 
sum_flux!(du, o::Organism) = begin
    offset_apply!(sum_flux!, du, o.organs, 0)
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
keep_records!(o::Organism, t) = begin
    apply(keep_records!, o.organs, o.records, t)
    o
end
keep_records!(organ, records::Records, t) = 
    organ.vars, organ.J, organ.J1 = records.vars[t], records.J[t], records.J1[t]
keep_records!(organ, records::Void, t) = nothing

"Handle dual number or other types in du if needed"
check_du_type(du, o) = begin
    du_fill = oneunit_flux(o.organs[1].params, o.organs[1].state)
    if typeof(du_fill) == eltype(du) 
        du
    else
        fill(du_fill, 12) 
    end
end
