" Surface area scaling rules "
abstract type AbstractShape end

struct Isomorph <: AbstractShape end 

@columns struct V0morph{Mo} <: AbstractShape
    Vd::Mo  | 4.0 | mol | [0.0, 1000.0] | _ | "reference"
end

@columns struct V1morph{Mo} <: AbstractShape
    Vd::Mo  | 4.0 | mol | [0.0, 1000.0] | _ | "reference"
end

@columns struct V1V0morph{Mo,B} <: AbstractShape
    Vd::Mo   | 4.0 | mol | [0.0, 1000.0] | _ | "reference"
    Vmax::Mo | 4.0 | mol | [0.0, 1000.0] | _ | "reference"
    β::B     | 4.0 | mol | [0.0, 10.0]   | _ | "reference"
end

@columns struct Plantmorph{Mo} <: AbstractShape
    M_Vref::Mo     | 0.1 | mol | [0.00001, 20.0]  | true | "Scaling reference"
    M_Vscaling::Mo | 1.0 | mol | [0.0001, 2000.0] | true | "Scaling mass"
end

update_shape!(o, u) = set_shape!(o, shape_correction(shape_pars(o), u.V))

shape_correction(f::Nothing, V) = 1
shape_correction(f::Isomorph, V) = 1
shape_correction(f::V0morph, V) = (V / f.Vd)^(-2//3)
shape_correction(f::V1morph, V) = (V / f.Vd)^(1//3)
shape_correction(f::V1V0morph, V) = (V / f.Vd)^(1//3 - (V/f.Vmax)^f.β)
shape_correction(f::Plantmorph, V) = (V / f.M_Vref)^(-V / f.M_Vscaling)

