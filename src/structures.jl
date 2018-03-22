using Base.tail

export scale_time_dependent_params!get_num_states, get_state_names, get_total_states

get_state_names(u::FieldVector) = fieldnames(u)
get_num_states(p::S) where S<:AbstractStructuredSettings = 
    length(fieldnames(p.structures[1].u))
get_state_names(p::S) where S<:AbstractStructuredSettings = 
    get_state_names(p.structures[1].u)
get_total_states(p::S) where S<:AbstractSettings = length(p.u0)


# -----------------------------------------------------------------------------
# Outer constructor for DEBStructure
DEBStructure(name::Symbol, param_specs::ParamSpecs, params, flags::DEBFlags,
             functions, settings) = begin
    (tmin, tmax) = ceil.(Int, settings[:tspan])
    timerange = tmin:tmax
    arraylength = tmax - tmin + 1

    rates = fill(0.0, arraylength)
    param_ids = [p.id for (k, p) in param_specs]

    A = 0.0
    init_params = deepcopy(params)
    u = init_state(settings[:state_type])
    num_structures = floor(Int64, length(settings[:u0])/length(u))
    assim_state = PhotoState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                             0.0, 0.0, 0, 0, false)

    state_names = fieldnames(u)
    state1_names = get_state1_names(u)
    Jbase = build_axis(state_names, TRANS, arraylength, timerange)
    J1base = build_axis(state1_names, TRANS1, arraylength, timerange)
    J = split_flux(Jbase, 1)
    J1 = split_flux(J1base, 1)
    structure = DEBStructure(name, param_specs, param_ids, init_params, params, flags, 
                             functions, u, A, Jbase, J1base, J, J1, rates, assim_state)
    return structure
end


function split_flux(base::AbstractArray, t::Int)
    view(base, :, :, t)
end

function build_axis(x, y, len, timerange)
    AxisArray(zeros(Float64, length(x), length(y), len),
              Axis{:state}(x),
              Axis{:transformations}(y),
              Axis{:time}(timerange))
end


# -----------------------------------------------------------------------------
# Structure update functions. 
# Recursion is used instead of loops for type stability.

split_state!(ss::Tuple{<:DEBStructure,Vararg}, offset::Int, u::A where A)::Void = begin
    states = length(fieldnames(ss[1].u))
    for i = 1:states
        ss[1].u[i] = u[i + offset]
    end
    split_state!(tail(ss), offset + states, u)
end
split_state!(::Tuple{}, ::Int, _)::Void = nothing


initialise_params!(s::S where S) = s.params .= s.init_params

scale_time_dependent_params!(s::S where S, timestep_days::Float64) = begin 
    for id in s.flags.time
        # TODO this is dangerous, it should use the param_specs value.
        s.init_params[id] *= timestep_days
    end
end

set_current_flux!(ss::Tuple{<:DEBStructure,Vararg}, t::Int)::Void = begin
  s = ss[1]
  s.J = split_flux(s.Jbase, t + 1)
  s.J1 = split_flux(s.J1base, t + 1)

  set_current_flux!(tail(ss), t)
end
set_current_flux!(::Tuple{}, ::Int)::Void = nothing


sum_flux!(du::AbstractArray, p::P where P<:AbstractStructuredSettings)::Void = begin
    ss = p.structures
    num_structures = length(p.structures)
    trans = length(TRANS)
    states = length(fieldnames(ss[1].u))
    sum_flux!(du, ss, num_structures, 0, states, trans)
    return nothing
end
function sum_flux!(du::AbstractArray, ss::Tuple{<:DEBStructure,Vararg},
                   num_structures::Int, offset::Int, states::Int, trans::Int)::Void
    J = ss[1].J
    for x = 1:states
        sum = 0.0
        for y = 1:trans
            sum += J[x, y]
        end
        du[x + offset] = sum
    end
    sum_flux!(du, tail(ss), num_structures, offset + states, states, trans)
end
sum_flux!(::AbstractArray, ::Tuple{}, ::Int, ::Int, ::Int, ::Int)::Void = nothing
