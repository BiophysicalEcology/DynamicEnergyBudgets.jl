"""
Maturity formulations allocate a fraction of 
resources to maturity and reproduction.
"""
abstract type AbstractMaturity end

"""
    maturity!(o, u)
    maturity!(f, o, u)

Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
function maturity! end
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Nothing, o, u) = nothing

"""
    Maturity(j_E_mat_mai, κmat, threshold)

A maturity model seperated to make maturity modeling optional.

$(FIELDDOCTABLE)
"""
@columns struct Maturity{MoMoD,F,Mo} <: AbstractMaturity
    # Field            | Default         | Unit            | Bounds       | Log  | Description
    j_E_mat_mai::MoMoD | 0.001           | mol*mol^-1*d^-1 | (0.0, 0.1)   | _    | "Spec maturity maint costs "
    κmat::F            | 0.05            | _               | (0.0, 1.0)   | _    | "Reserve flux allocated to development/reprod."
    threshold::Mo      | 1.0             | mol             | (1e-3, 20.0) | true | "Structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  
end

# n_N_M::MoMo      | 0.05            | mol*mol^-1      | [0.0, 1.0]   | _    | "N/C use for maturity"
# w_M::GMo         | 25.0            | g*mol^-1        | [0.0, 1.0]   | _    | "Mol-weight of shoot maturity reserve:"

maturity!(f::Maturity, o, u) = begin
    seedset = u[:M] > f.threshold ? u[:M] : zero(u[:M]) 
    mat = κmat(f) * E_ctb(o)
    flux(o)[:M,:gro] = mat - seedset * unit(mat)/unit(seedset)
    mat_mai = f.j_E_mat_mai * tempcorrection(o) * u[:V]
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
end

"""
    κmat(maturity_pars::Maturity)

Kappa parameter for maturity.
"""
κmat(maturity_pars::Maturity) = maturity_pars.κmat
