"""
Outer construtor for defaults, without J and J1
"""
Organ(; state = StatePVMCNE(),
        name = :Shoot,
        params = Params(),
        shared = SharedParams(),
        vars = Vars(),
        time = 0u"hr":1u"hr":1000u"hr"
     ) = begin

    J = build_flux(fieldnames(state), TRANS)
    J1 = build_flux(get_state1_names(state), TRANS1)
    Organ(state, name, params, shared, vars, J, J1, time, varsrecord, Jrecord, J1record)
end

"""
Outer construtor for defaults, without J and J1
"""
Organism(; records = nothing,
           shared = SharedParams(),
           nodes = (Organ(time=time), 
                    Organ(time=time,
                          params=Params(assimilation=Kooijman_NH4_NO3Assimilation()),
                          vars=Vars(assimilation=NitrogenVars())))) = begin
    records = []
    for organ in nodes
        organ.shared = shared
        push!(records, Records(organ, time))
    end
    Organism(tuple(records...), shared, nodes)
end

" build all required records with the length of the current timespan "
Records(o::Organ, time) = begin
    Jrec = build_record(build_J(o.state), time)
    J1rec = build_record(build_J1(o.state), time)
    varsrec = build_record(vars, time)
    Records(varsrec, Jrec, J1rec)
end

build_J(state) = build_flux(fieldnames(state), TRANS)
build_J1(state) = build_flux(get_state1_names(state), TRANS1)

build_flux(x, y) = begin
    a = fill(0.0u"mol*hr^-1", length(x), length(y))
    AxisArray(a, Axis{:state}(x), Axis{:transformations}(y))
end

build_record(a, time) = AxisArray([deepcopy(a) for i = 1:length(time)], Axis{:time}(time))

" update vars and flux to record for time t "
set_cur_records!(organ, records::Records, t) = begin
    organ.vars = records.vars[t] 
    organ.J = records.J[t]
    organ.J1 = records.J1[t]
end
" Void records so we just stay on the same record "
set_cur_records!(organ, records::Void, t) = nothing

" copy the diffeq state to organs "
setstate!(organ, u, offset::Int) = begin
    for i in 1:length(organ.state) 
        organ.state[i] = u[i+offset]
    end
    offset + length(organ.state)
end

" sum flux matrix "
sumflux!(du, organ, offset::Int) = begin
    for i in 1:size(organ.J, 1) 
        du[i+offset] = sum(organ.J[i,:])
    end
    offset + length(organ.state)
end

