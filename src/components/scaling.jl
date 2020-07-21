"""
Surface area scaling rules
"""
abstract type AbstractScaling end

"""
    scaling_correction(f::AbstractScaling, V)

Calculate the shape/scaling correction coefficient 
from the current mass `V`.
"""
function scaling_correction end
scaling_correction(f::Nothing, V) = 1

"""
    Isomorph()
"""
struct Isomorph <: AbstractScaling end

scaling_correction(f::Isomorph, V) = 1

"""
    V0morph(Vd)

$(FIELDDOCTABLE)
"""
@columns struct V0morph{Mo} <: AbstractScaling
#   Field    | Def | Unit | Bounds        | Log | Description
    Vd::Mo   | 4.0 | mol  | (0.0, 1000.0) | _   | "reference"
end

scaling_correction(f::V0morph, V) = (V / f.Vd)^(-2//3)

"""
    V1morph(Vd)

$(FIELDDOCTABLE)
"""
@columns struct V1morph{Mo} <: AbstractScaling
#   Field    | Def | Unit | Bounds        | Log | Description
    Vd::Mo   | 4.0 | mol  | (0.0, 1000.0) | _   | "reference"
end

scaling_correction(f::V1morph, V) = (V / f.Vd)^(1//3)

"""
    V1V0morph(Vd, Vmax, β)

$(FIELDDOCTABLE)
"""
@columns struct V1V0morph{Mo,B} <: AbstractScaling
#   Field    | Def | Unit | Bounds        | Log | Description
    Vd::Mo   | 4.0 | mol  | (0.0, 1000.0) | _   | "reference"
    Vmax::Mo | 4.0 | mol  | (0.0, 1000.0) | _   | "reference"
    β::B     | 4.0 | mol  | (0.0, 10.0)   | _   | "reference"
end

scaling_correction(f::V1V0morph, V) = (V / f.Vd)^(1//3 - (V/f.Vmax)^f.β)

"""
    Plantmorph(M_Vref, M_Vscaling)

Plant morph formulation from DEBtool.

$(FIELDDOCTABLE)
"""
@columns struct Plantmorph{Mo} <: AbstractScaling
#   Field          | Def | Unit | Bounds           | Log  | Description
    M_Vref::Mo     | 0.1 | mol  | (0.00001, 20.0)  | true | "Scaling reference"
    M_Vscaling::Mo | 1.0 | mol  | (0.0001, 2000.0) | true | "Scaling mass"
end

scaling_correction(f::Plantmorph, V) = (V / f.M_Vref)^(-V / f.M_Vscaling)

update_scaling!(o, u) = set_scaling!(o, scaling_correction(scaling_pars(o), u[:V]))
