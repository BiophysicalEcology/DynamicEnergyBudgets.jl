abstract type AbstractDEBCore end

"""
    DEBCore(y_V_E, y_E_C, y_E_N, n_N_V, n_N_E, w_V)

Core DEB model parameters.

$(FIELDDOCTABLE)
"""
@columns struct DEBCore{MoMoD,MoMo,GMo} <: AbstractDEBCore
    # Field        | Default | Unit            | Bounds       | Log   | Description
    j_E_mai::MoMoD | 0.01    | mol*mol^-1*d^-1 | (1e-4, 1.0)  | true  | "Specific somatic maintenance costs"
    y_V_E::MoMo    | 0.7     | _               | (0.0, 1.0)   | _     | "Yield from reserve to structure"
    y_E_C::MoMo    | 0.7     | _               | (1e-6, 1.0)  | false | "Yield from C-reserve to general reserve"
    y_E_N::MoMo    | 30.0    | _               | (1.0, 50.0)  | false | "Yield from N-reserve to general reserve"
    n_N_V::MoMo    | 0.03    | _               | (0.0, 0.1)   | _     | "Nitrogen per Carbon in structure"
    n_N_E::MoMo    | 0.025   | _               | (0.0, 0.1)   | _     | "Nitrogen per Carbon in reserve"
    w_V::GMo       | 25.0    | g*mol^-1        | (15.0, 40.0) | _     | "Mol-weight of shoot structure"
    # w_N::GMo     | 25.0    | g*mol^-1        | (15.0, 40.0) | _     | "Mol-weight of shoot N-reserve"
    # w_C::GMo     | 25.0    | g*mol^-1        | (12.0, 40.0) | _     | "Mol-weight of shoot C-reserve"
    # w_E::GMo     | 25.0    | g*mol^-1        | (15.0, 40.0) | _     | "Mol-weight of shoot reserve"
end

for fn in fieldnames(DEBCore)
    @eval $fn(p::DEBCore) = p.$fn
end

"""
    growth!(o::AbstractOrgan, u)
    growth!(p::DEBCore, o::AbstractOrgan, u)

Allocates reserves to growth flux, generalised for any number of reserves.

Where `o` is the `Organ`, and `u` is the current state parameters
"""
function growth! end
growth!(o, u) = growth!(core_pars(o), o, u)
growth!(p::DEBCore, o, u) = begin
    J = flux(o) 
    J[:V,:gro] = growth = rate(o) * u[:V]
    drain = (1 / y_V_E(p)) * growth 
    product = growth_production!(o, growth)
    reserve_drain!(o, :gro, drain)
end

"""
    maintenence!(o::AbstractOrgan, u)
    maintenence!(p::DEBCore, o::AbstractOrgan, u)

Allocates reserve drain due to maintenance, generalised for any number of reserves.

Maintenance is temperature dependent.

Where `o` is the `Organ`, and `u` is the current state parameters
"""
function maintenance! end
maintenence!(o, u) = maintenence!(core_pars(o), o, u)
maintenence!(p::DEBCore, o, u) = begin
    drain = j_E_mai(p) * tempcorrection(o) * u[:V]
    reserve_drain!(o, :mai, drain)
    maintenance_production!(o, u)
end

corrected_j_E_mai(o) = j_E_mai(o) * tempcorrection(o)
