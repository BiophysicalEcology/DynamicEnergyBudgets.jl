# Translocation of reserve rejected from the SU during catabolism

abstract type AbstractRejection end

"""
    translocate_rejected!(source, dest, prop)

Reallocate state rejected from synthesizing units.

TODO: add a 1-organs method. How does this interact with assimilation?
"""
translocate_rejected!(source, dest, prop) =
    translocate_rejected!(rejection_pars(source), source, dest, prop)
translocate_rejected!(f::Nothing, source, dest, prop) = nothing

"""
    DissipativeRejection(y_EC_ECT, y_EN_ENT)

Substrate rejected from synthesizing units during catabolism is returned to
reserve, but with some fraction of loss specified by yield parameters.
"""
@columns struct DissipativeRejection{MoMo} <: AbstractRejection
#   Field          | Default | Unit       | Bounds     | Log | Description
    y_EC_ECT::MoMo | 1.0     | mol*mol^-1 | (0.0, 1.0) | _   | "yield of translocated C-reserve"
    y_EN_ENT::MoMo | 1.0     | mol*mol^-1 | (0.0, 1.0) | _   | "yield of Translocated N-reserve"
end

translocate_rejected!(f::DissipativeRejection, source, dest, prop) = begin
    Js, Jd = flux(source), flux(dest)
    Jd[:C,:tra] = -Js[:C,:rej] * f.y_EC_ECT
    Jd[:N,:tra] = -Js[:N,:rej] * f.y_EN_ENT
    nothing
end

"""
    LosslessRejection)

Parameterless rejection where substrate rejected from synthesizing units
during catabolism is returned to reserve without loss.
"""
struct LosslessRejection <: AbstractRejection end

translocate_rejected!(rejected::LosslessRejection, source, dest, prop) = begin
    Js, Jd = flux(source), flux(dest)
    Jd[:C,:tra] = -Js[:C,:rej]
    Jd[:N,:tra] = -Js[:N,:rej]
    nothing
end


# Active translocation. Not required in practice, as
# translocation of rejected reserves gives resonable behaviour.

abstract type AbstractTranslocation end

@mix @flattenable @columns struct Trans{F}
#   Field          | Flat | Default  | Unit       | Bounds    | Log | Description
    κtra::F        | true | 0.6      | _          | (0.0,1.0) | _   | "Reserve flux allocated to translocation"
end

κtra(trans_pars::AbstractTranslocation) = trans_pars.κtra

@mix @flattenable @columns struct MultiTrans{D,P}
    destnames::D   | false| (:leaf,) | _          | _         | _   | "The organ/s translocated to"
    proportions::P | true | (1.0,)   | _          | (0.0,1.0) | _   | "The proportion of translocation sent in the first translocation. Only for inetermediaries. nothing = 100%"
end

@mix @flattenable @columns struct DissTrans{MoMo}
    y_E_ET::MoMo   | true | 0.8      | mol*mol^-1 | (0.0,1.0) | _   | "yeild of translocated reserve:"
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


abstract type AbstractDissipativeTranslocation <: AbstractTranslocation end

"""
    DissipativeTranslocation(κtra, y_E_ET)

Translocation with dissipative losses to the environment.
"""
@Trans @DissTrans struct DissipativeTranslocation{} <: AbstractDissipativeTranslocation  end

translocate!(p::DissipativeTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, :tra, trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb]
    o1.J[:E,:tra] += transx * y_E_ET(o2)

    # loss = transx * (1 - y_E_ET(o2))
    # reserve_loss!(o2, loss)
    # conversion_loss!(o2, transx * y_E_ET(o2), n_N_E(o2))
    nothing
end


abstract type AbstractLosslessTranslocation <: AbstractTranslocation end

"""
    DissipativeTranslocation(κtra, y_E_ET)

Perfect translocation between structures.
"""
@Trans struct LosslessTranslocation{} <: AbstractLosslessTranslocation end

translocate!(p::LosslessTranslocation, o1, o2, prop) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb]
    reserve_drain!(o1, :tra, trans)

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
