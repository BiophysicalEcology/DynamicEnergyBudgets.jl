abstract type AbstractMaintenance end

@columns struct Maintenance{MoMoD} <: AbstractMaintenance
    # Field        | Def  | Unit            | Bounds      | Log  | Description
    j_E_mai::MoMoD | 0.01 | mol*mol^-1*d^-1 | (1e-4, 1.0) | true | "Specific somatic maintenance costs"
end

"""
Allocates reserve drain due to maintenance.
"""
maintenence!(o, u) = maintenence!(maintenance_pars(o), o, u)
maintenence!(f::Maintenance, o, u) = begin
    drain = j_E_mai(f) * tempcorrection(o) * u.V
    maint_prod = maintenance_production!(o, u)
    reserve_drain!(o, Val(:mai), drain)
    # reserve_loss!(o, drain - maint_prod) # all maintenance is loss
end

j_E_mai(f::Maintenance) = f.j_E_mai
