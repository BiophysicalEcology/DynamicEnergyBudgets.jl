" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns struct SqrtAllometry{M} <: AbstractAllometry
    # Field   | Default | Unit | Pror            | Limits      | Log  | Description
    scalar::M | 0.1     | m    | Gamma(2.0, 0.2) | [1e-2, 1.0] | true | "Scalar"
end

@columns struct Allometry{M} <: AbstractAllometry
    # Field | Default | Unit | Prior           | Limits         | Log  | Description
    β::B    | 0.1     | m    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Scalar"
    α::A    | 0.1     | _    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Exponent"
end

@columns struct FixedAllometry{M} <: AbstractAllometry
    # Field   | Default | Unit | Prior           | Limits         | Log  | Description
    height::B | 1.0     | m    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Fixed height or depth"
end


update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(a::Nothing, o, u) = nothing
update_height!(a, o, u) = set_height!(o, allometric_height(a, o, u))

allometric_height(f::SqrtAllometry, o, u) = begin
    units = unit(u.V * w_V(o))
    sqrt((u.V * w_V(o)) / units) * f.scalar
end

allometric_height(f::Allometry, o, u) = begin
    units = unit(u.V * w_V(o))
    f.β * (u.V * w_V(o)/units)^f.α
end

allometric_height(f::Allometry, o, u) = f.height
