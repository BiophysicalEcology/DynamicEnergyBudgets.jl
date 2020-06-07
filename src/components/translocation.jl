abstract type AbstractTranslocation end
abstract type AbstractDissipativeTranslocation <: AbstractTranslocation end
abstract type AbstractLosslessTranslocation <: AbstractTranslocation end

@mix @flattenable @columns struct Trans{F}
#   Field          | Flat  | Default  | Unit | Bounds |
    κtra::F        | true  | 0.6      | _ | [0.0,1.0] | _ | "Reserve flux allocated to translocation"
end

@mix @flattenable @columns struct MultiTrans{D,P}
    destnames::D   | false | (:leaf,) | _ | _         | _ | "The organ/s translocated to"
    proportions::P | true  | (1.0,)   | _ | [0.0,1.0] | _ | "The proportion of translocation sent in the first translocation. Only for inetermediaries. nothing = 100%"
end

@mix @flattenable @columns struct DissTrans{MoMo}
    y_E_ET::MoMo   | true | 0.8       | mol*mol^-1   | [0.0,1.0] | _ | "Translocated reserve:"
end



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



@Trans @DissTrans struct DissipativeTranslocation{} <: AbstractDissipativeTranslocation  end

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


@Trans struct LosslessTranslocation{} <: AbstractLosslessTranslocation end

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


# Not implemented yet
@Trans @MultiTrans struct LosslessMultipleTranslocation{} <: AbstractLosslessTranslocation end
@Trans @MultiTrans @DissTrans struct DissipativeMultipleTranslocation{} <: AbstractDissipativeTranslocation end

"""
Translocation occurs between adjacent organs.
This function is identical both directiono, so on represents
whichever is not the current organs.

Will not run with less than 2 organs.
"""
translocation!(organs::Tuple{T1,T2}) where {T1, T2} = begin
    translocate_rejected!(organs[1], organs[2], 1.0)
    translocate_rejected!(organs[2], organs[1], 1.0)
    translocate!(organs[1], organs[2], 1.0)
    translocate!(organs[2], organs[1], 1.0)
end
# translocation!(organs::Tuple) = translocation!(organs...)
translocation!(organs::Tuple{T}) where T = nothing
translocation!(organs::Tuple{}) = nothing



abstract type AbstractRejection end

"""
Reallocate state rejected from synthesizing units.
TODO: add a 1-organs method. How does this interact with assimilation?  
"""
translocate_rejected!(source, dest, prop) = translocate_rejected!(rejection_pars(source), source, dest, prop)
translocate_rejected!(f::Nothing, source, dest, prop) = nothing

"""
    DissipativeRejection(y_EC_ECT, y_EN_ENT)
"""
@columns struct DissipativeRejection{MoMo} <: AbstractRejection
    y_EC_ECT::MoMo       | 1.0             | mol*mol^-1      | [0.0,1.0] | _ | "Translocated C-reserve"
    y_EN_ENT::MoMo       | 1.0             | mol*mol^-1      | [0.0,1.0] | _ | "Translocated N-reserve"
end

translocate_rejected!(f::DissipativeRejection, source, dest, prop) = begin
    Js, J1s, Jd = flux(source), flux1(source), flux(dest)
    transC = J1s[:C,:rej] # * (1 - κEC(o))
    transN = J1s[:N,:rej] # * (1 - κEN(o))
    Js[:C,:rej] = -transC
    Js[:N,:rej] = -transN
    Jd[:C,:tra] = f.y_EC_ECT * transC
    Jd[:N,:tra] = f.y_EN_ENT * transN
    # J1[:C,:los] += transC * (1 - y_EC_ECT(o)) + transN * (1 - y_EN_ENT(o))
    # J1[:N,:los] += (transC * (1 - y_EC_ECT(o)), transN * (1 - y_EN_ENT(o))) ⋅ (n_N_EC(o), n_N_EN(o))
    nothing
end

struct LosslessRejection <: AbstractRejection end

translocate_rejected!(rejected::LosslessRejection, source, dest, prop) = begin
    Js, J1s, Jd = flux(source), flux1(source), flux(dest)
    transC = J1s[:C,:rej] # * (1 - κEC(o))
    transN = J1s[:N,:rej] # * (1 - κEN(o))
    Js[:C,:rej] = -transC
    Js[:N,:rej] = -transN
    Jd[:C,:tra] = transC
    Jd[:N,:tra] = transN
    nothing
end


κtra(trans_pars::AbstractTranslocation) = trans_pars.κtra
