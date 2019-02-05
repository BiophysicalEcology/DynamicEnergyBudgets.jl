abstract type AbstractTranslocation end
abstract type AbstractDissipativeTranslocation <: AbstractTranslocation end
abstract type AbstractLosslessTranslocation <: AbstractTranslocation end

@mix @columns @flattenable struct Trans{F}
    κtra::F        | 0.6     | _   | Beta(2.0, 2.0)  | [0.0,1.0]      | _ | "Reserve flux allocated to translocation"
end

@mix @columns @flattenable struct MultiTrans{D,P}
    destnames::D   | false | (:leaf,) | _ | _              | _         | _ | "The organ/s translocated to"
    proportions::P | true  | (1.0,)   | _ | Beta(2.0, 2.0) | [0.0,1.0] | _ | "The proportion of translocation sent in the first translocation. Only for inetermediaries. nothing = 100%"
end

@mix @columns @flattenable struct DissTrans{MoMo}
    y_E_ET::MoMo   | true | 0.8       | mol*mol^-1      | Beta(2.0, 2.0)  | [0.0,1.0]    | _ | "Translocated reserve:"
end

# Dissipation types
@Trans @MultiTrans @DissTrans struct DissipativeMultipleTranslocation{} <: AbstractDissipativeTranslocation end
@Trans @DissTrans struct DissipativeTranslocation{} <: AbstractDissipativeTranslocation  end
@Trans @MultiTrans struct LosslessMultipleTranslocation{} <: AbstractLosslessTranslocation end
@Trans struct LosslessTranslocation{} <: AbstractLosslessTranslocation end

# Recurse through all organs. A loop would not be type-stable.
# translocation!(organs::Tuple, destorgans::Tuple) = begin
#     props = buildprops(organs[1])
#     translocation!(organs[1], destorgans, organs[1].params.translocation.destnames, props)
#     translocation!(tail(organs), destorgans)
# end
# translocation!(organs::Tuple{}, destorgans::Tuple) = nothing
# translocation!(organ::Organ, destorgans::Tuple, destnames::Symbol, props) =
#     translocation!(organ, destorgans, (destnames,), props)
# # Translocate to organs with names in the destnames list
# translocation!(organ::Organ, destorgans::Tuple, destnames, props) = begin
#     for i = 1:length(destnames)
#         if destorgans[1].params.name == destnames[i]
#             reuse_rejected!(organ, destorgans[1], props[i])
#             translocate!(organ, destorgans[1], props[i])
#             break
#         end
#     end
#     translocation!(organ, tail(destorgans), destnames, props)
# end
# translocation!(organ::Organ, destorgans::Tuple{}, destnames, props) = nothing

# Add the last remainder proportion (so that its not a model parameter)
# buildprops(o::Organ) = buildprops(o.params.translocation.proportions)
# buildprops(x::Nothing) = (1.0)
# buildprops(x::Number) = (x, 1 - x)
# buildprops(xs::Tuple) = (xs..., 1 - sum(xs))

"""
Versions for E, CN and CNE reserves.

Translocation is occurs between adjacent organs.
This function is identical both directiono, and ox represents
whichever is not the current organs. Will not run with less than 2 organs.
"""
translocate!(o1, o2, prop) = translocate!(trans_pars(o1), o1, o2, prop)
translocate!(p::Nothing, o1, o2, prop) = nothing
translocate!(p::AbstractDissipativeTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx * y_E_ET(o2)

    # loss = transx * (1 - y_E_ET(o2))
    # reserve_loss!(o2, loss)
    # conversion_loss!(o2, transx * y_E_ET(o2), n_N_E(o2))
    nothing
end

translocate!(p::AbstractLosslessTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, Val(:tra), trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx

    # conversion_loss!(o2, transx, n_N_E(o2))
    nothing
end
