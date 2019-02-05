"""
Generalised reserve loss to track carbon.

Loss is distributed between general and C and N reserves by the fraction θE
"""
@inline reserve_loss!(o, loss) = nothing #reserve_loss!(o.params, o, loss)
@inline reserve_loss!(::HasCNE, o, loss) = begin
    ee = loss * θ # fraction of loss from E reserve
    ecn = loss - ee # fraction on loss from C and N reserves
    ec = ecn/y_E_EC(o)
    en = ecn/y_E_EN(o)
    o.J1[:C,:los] += ec + en + ee
    o.J1[:N,:los] += (ec, en, ee) ⋅ (n_N_EC(o), n_N_EN(o), n_N_E(o))
    nothing
end
@inline reserve_loss!(::HasCN, o, loss) = begin
    ec = loss/y_E_EC(o)
    en = loss/y_E_EN(o)
    o.J1[:C,:los] += ec + en
    o.J1[:N,:los] += (ec, en) ⋅ (n_N_EC(o), n_N_EN(o))
    nothing
end

@inline conversion_loss!(o, loss, dest_n_N) = nothing #conversion_loss!(o.params, o, loss, dest_n_N)
@inline conversion_loss!(::HasCNE, o, loss, dest_n_N) = begin
    ecn = loss * (1 - θE(o)) # fraction on loss from C and N reserves
    ec = ecn/y_E_EC(o)
    en = ecn/y_E_EN(o)
    o.J1[:C,:los] += ec + loss * (θE(o) - 1) # + en
    o.J1[:N,:los] += (ec, en, loss * (θ - dest_n_N/n_N_E(o))) ⋅ (n_N_EC(o), n_N_EN(o), n_N_E(o))
end
@inline conversion_loss!(::HasCN, o, loss, dest_n_N) = begin
    ec = loss/y_E_EC(o)
    en = loss/y_E_EN(o)
    o.J1[:C,:los] += ec + loss # + en
    o.J1[:N,:los] += (ec, en) ⋅ (n_N_EC(o), n_N_EN(o))
end

stoich_merge_losses(Jc1, Jn1, Jc2, Jn2, JEcn, n_c, n_n, n_Ecn) = begin
    lossa = Jc1 - Jc2 + Jn1 - Jn2 - JEcn
    lossb = (Jc1 - Jc2, Jn1 - Jn2, -JEcn) ⋅ (n_c, n_n, n_Ecn)
    lossa, lossb
end
