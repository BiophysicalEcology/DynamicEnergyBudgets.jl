# Outer construtor for defaults, without J and J1
Organ(; state = StatePVMCNE(),
        name = :Shoot,
        params = Params(),
        shared = SharedParams(),
        vars = Vars(),
        time = 0u"hr":1u"hr":1000u"hr"
     ) = begin

    varsrecord, Jrecord, J1record = timespan_records(time, state, vars)
    J, J1 = Jrecord[1], J1record[1]
    Organ(state, name, params, shared, vars, J, J1, time, varsrecord, Jrecord, J1record)
end


Organism(; time = 0u"hr":1u"hr":1000u"hr",
           shared = SharedParams(),
           nodes = (Organ(time=time), 
                    Organ(time=time,
                          params=Params(assimilation=Kooijman_NH4_NO3Assimilation()),
                          vars=Vars(assimilation=NitrogenVars())))) = begin
    for organ in nodes
        organ.shared = shared
        set_timespan(organ, time)
    end
    Organism(time, shared, nodes)
end

build_axis(x, y, time) = begin
    a = fill(0.0u"mol*hr^-1", length(x), length(y))
    axis = AxisArray(a, Axis{:state}(x), Axis{:transformations}(y))
    AxisArray([deepcopy(axis) for i = 1:length(time)], Axis{:time}(time))
end


timespan_records(time, state, vars) = begin 
    Jrecord = build_axis(fieldnames(state), TRANS, time)
    J1record = build_axis(get_state1_names(state), TRANS1, time)
    a = [deepcopy(vars) for t in time]
    varsrecord = AxisArray(a, Axis{:time}(time))
    varsrecord, Jrecord, J1record 
end

set_timespan(o::Organ, time) = begin
    o.varsrecord, o.Jrecord, o.J1record = timespan_records(time, o.state, o.vars)
    o.J = o.Jrecord[1]
    o.J1 = o.J1record[1]
end
set_timespan(o::Organism, time) = set_timespan(o.nodes, time)
set_timespan(o::Scenario, time) = set_timespan(o.nodes, time)
