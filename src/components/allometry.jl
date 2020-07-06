" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(f::Nothing, o, u) = nothing
update_height!(f::AbstractAllometry, o, u) =
    set_height!(o, allometric_height(f, mass(o, u)))

"""
    SqrtAllometry(β0, β)

Height is given by the square root of mass above the initial mass `β0`,
multiplied by the scalar `β1`.
"""
@flattenable @columns struct SqrtAllometry{B0,B1} <: AbstractAllometry
    # Field       | Flatn | Default | Unit | Bounds         | Log  | Description
    β0::B0        | false | 1e-4*24 | g    | (1e-6, 10.00)  | true | "Intercept. Mass at Om"
    β1::B1        | true  | 0.1     | m    | (1e-2, 1.0)    | true | "Scalar for conversion to meters"
end

allometric_height(f::SqrtAllometry, mass) =
    sqrt((mass - f.β0) / unit(f.β0)) * f.β1

"""
    Allometry(β0, β, α)

Simple allometric relationship between mass and height
"""
@flattenable @columns struct Allometry{B0,B1,A} <: AbstractAllometry
    β0::B0        | false | 1e-4*24 | g    | (1e-6, 10.00)  | true | "Intercept. Mass at Om"
    β1::B1        | true  | 0.1     | m    | (1e-5, 10.00)  | true | "Scalar for conversion to meters"
    α::A          | true  | 0.1     | _    | (1e-3, 100.00) | true | "Exponent relating mass to vertical dimension"
end

allometric_height(f::Allometry, mass) =
    f.β1 * ((max(zero(mass), mass - f.β0)) / unit(f.β0))^f.α

"""
    FixedAllometry(height)

Height is fixed at `height`, independent of mass.
"""
@flattenable @columns struct FixedAllometry{M} <: AbstractAllometry
    height::M     | true  | 1.0     | m    | (1e-3, 100.00) | true | "Fixed height or depth"
end

allometric_height(f::FixedAllometry, mass) = f.height

mass(o, u) = u[:V] * w_V(o)
