" Allometry. Scaling rules to relate size to mass. "
abstract type AbstractAllometry end

@columns struct SqrtAllometry{M} <: AbstractAllometry
    # Field         | Default  | Unit | Prior           | Limits        | Log | Description
    allometry::M    | 0.1      | m    | Gamma(2.0, 0.2) | [0.0, 0.20]   | _   | "Allometric height/depth scaling"
end
