" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns @flattenable struct SqrtAllometry{M,I} <: AbstractAllometry
    # Field       | Flatn | Default | Unit | Pror            | Limits      | Log  | Description
    scalar::M     | true  | 0.1     | m    | Gamma(2.0, 0.2) | [1e-2, 1.0]    | true  | "Scalar"
    inputunits::I | false | g       | _    | _               | _              | false | "Input units"
end

@columns @flattenable struct Allometry{B,A,I} <: AbstractAllometry
    β::B          | true  | 0.1     | m    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true  | "Scalar"
    α::A          | true  | 0.1     | _    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true  | "Exponent"
    inputunits::I | false | g       | _    | _               | _              | false | "Input units"
end

@columns struct FixedAllometry{M} <: AbstractAllometry
    height::M             | 1.0     | m    | Gamma(2.0, 0.2) | [1e-3, 100.00] | true | "Fixed height or depth"
end


update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(f::Nothing, o, u) = nothing
update_height!(f, o, u) = set_height!(o, allometric_height(f, mass(o, u)))

allometric_height(f::SqrtAllometry, mass) = sqrt(mass/f.inputunits) * f.scalar
allometric_height(f::Allometry, mass) = f.β * (mass/f.inputunits)^f.α
allometric_height(f::FixedAllometry, mass) = f.height
