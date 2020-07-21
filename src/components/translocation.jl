# Translocation of reserve rejected from the SU during catabolism

abstract type PassiveTranslocation end

"""
    passive_translocation!(source, dest)

Realocate state rejected from synthesizing units.
"""
function passive_translocation! end
passive_translocation!(source, dest) =
    passive_translocation!(passivetrans_pars(source), source, dest)
passive_translocation!(f::Nothing, source, dest) = nothing

"""
    DissipativeRejection(y_EC_ECT, y_EN_ENT)

Substrate rejected from synthesizing units during catabolism is returned to
reserve, but with some fraction of loss specified by yield parameters.

$(FIELDDOCTABLE)
"""
@columns struct DissipativePassiveTranslocation{MoMo} <: PassiveTranslocation
#   Field          | Default | Unit       | Bounds     | Log | Description
    y_EC_ECT::MoMo | 1.0     | mol*mol^-1 | (0.0, 1.0) | _   | "yield of translocated C-reserve"
    y_EN_ENT::MoMo | 1.0     | mol*mol^-1 | (0.0, 1.0) | _   | "yield of Translocated N-reserve"
end

passive_translocation!(f::DissipativePassiveTranslocation, source, dest) = begin
    Js, Jd = flux(source), flux(dest)
    Jd[:C,:tra] += -Js[:C,:rej] * f.y_EC_ECT
    Jd[:N,:tra] += -Js[:N,:rej] * f.y_EN_ENT
    nothing
end

"""
    LosslessPassiveTranslocation()

Parameterless rejection where substrate rejected from synthesizing units
during catabolism is returned to reserve without loss.

$(FIELDDOCTABLE)
"""
struct LosslessPassiveTranslocation <: PassiveTranslocation end

passive_translocation!(rejected::LosslessPassiveTranslocation, source, dest) = begin
    Js, Jd = flux(source), flux(dest)
    Jd[:C,:tra] += -Js[:C,:rej]
    Jd[:N,:tra] += -Js[:N,:rej]
    nothing
end


# Active translocation. Not required in practice, as
# translocation of rejected reserves gives resonable behaviour.

abstract type ActiveTranslocation end

@mix @flattenable @columns struct Active{F}
#   Field          | Flat | Default  | Unit       | Bounds    | Log | Description
    κtra::F        | true | 0.6      | _          | (0.0,1.0) | _   | "Reserve flux allocated to translocation"
end

κtra(activetrans_pars::ActiveTranslocation) = activetrans_pars.κtra

"""
    active_translocation!(o1, o2) 

Translocation that actively moves a fraction of catabolised 
reserve between organs.
"""
function active_translocation! end
active_translocation!(o1, o2) = 
    active_translocation!(activetrans_pars(o1), o1, o2)
active_translocation!(p::Nothing, o1, o2) = nothing

"""
    DissipativeTranslocation(κtra, y_E_ET)

Translocation with dissipative losses to the environment.

$(FIELDDOCTABLE)
"""
@Active struct DissipativeActiveTranslocation{MoMo} <: ActiveTranslocation
    y_E_ET::MoMo   | true | 0.8      | mol*mol^-1 | (0.0,1.0) | _   | "yield of translocated reserve:"
end

active_translocation!(p::DissipativeActiveTranslocation, o1, o2) = begin
    # outgoing translocation
    trans = κtra(o1) * E_ctb(o1) / (1 - κsoma(o1))
    reserve_drain!(o1, :tra, trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb] / (1 - κsoma(o2))
    o1.J[:E,:tra] += transx * y_E_ET(o2)

    nothing
end


"""
    DissipativeTranslocation(κtra, y_E_ET)

Perfect translocation between structures.

$(FIELDDOCTABLE)
"""
@Active struct LosslessActiveTranslocation{} <: ActiveTranslocation end

active_translocation!(p::LosslessActiveTranslocation, o1, o2) = begin
    # outgoing translocation
    trans = κtra(o1) * o1.J1[:E,:ctb] / (1 - κsoma(o1))
    reserve_drain!(o1, :tra, trans)

    # incoming translocation
    transx = κtra(o2) * o2.J1[:E,:ctb] / (1 - κsoma(o1))

    o1.J[:E,:tra] += transx
    nothing
end

"""
    translocation!(organs::Tuple)

Translocation occurs between adjacent organs in 
both directions.

Both active and passive translocation are applied,
although `nothing` valued formulations do not translocate.
"""
function translocation! end
translocation!(organs::Tuple{T1,T2}) where {T1,T2} = begin
    passive_translocation!(organs[1], organs[2])
    passive_translocation!(organs[2], organs[1])
    active_translocation!(organs[1], organs[2])
    active_translocation!(organs[2], organs[1])
end
# translocation!(organs::Tuple) = translocation!(organs...)
translocation!(organs::Tuple{T}) where T = nothing
translocation!(organs::Tuple{}) = nothing
