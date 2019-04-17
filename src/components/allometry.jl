" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns @flattenable struct SqrtAllometry{B0,B1} <: AbstractAllometry
    # Field       | Flatn | Default | Unit | Pror            | Limits         | Log  | Description
    β0::B0        | false | 1e-4*24 | g    | Gamma(2.0, 0.2) | [1e-6, 10.00]  | true | "Intercept"
    β1::B1        | true  | 0.1     | m    | Gamma(2.0, 0.2) | [1e-2, 1.0]    | true | "Scalar"
end

@columns @flattenable struct Allometry{B0,B1,A} <: AbstractAllometry
    β0::B0        | false | 1e-4*24 | g    | Gamma(2.0, 0.2) | [1e-6, 10.00]  | true | "Intercept"
    β1::B1        | true  | 0.1     | m    | Gamma(2.0, 0.2) | [1e-5, 10.00]  | true | "Scalar"
    α::A          | true  | 0.1     | _    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Exponent"
end

@columns struct FixedAllometry{M} <: AbstractAllometry
    height::M             | 1.0     | m    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Fixed height or depth"
end


update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(f::Nothing, o, u) = nothing
update_height!(f, o, u) = set_height!(o, allometric_height(f, mass(o, u)))

allometric_height(f::FixedAllometry, mass) = f.height
allometric_height(f::SqrtAllometry, mass) = sqrt((mass - f.β0) / unit(f.β0)) * f.β1
allometric_height(f::Allometry, mass) = begin
    f.β1 * ((max(zero(mass), mass - f.β0))/unit(f.β0))^f.α
end
