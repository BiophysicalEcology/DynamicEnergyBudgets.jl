abstract type AbstractMaturity end

" Maturity parameters. Seperated to make maturity modeling optional, reducing complexity "
@columns struct ReproductionBuffer{MoMoD,F,Mo} <: AbstractMaturity
    # Field            | Default         | Unit            | Prior           | Limits      | Log  | Description
    # n_N_M::MoMo      | 0.05            | mol*mol^-1      | Gamma(2.0, 2.0) | [0.0, 1.0]  | _    | "N/C use for maturity"
    j_E_mat_mai::MoMoD | 0.001           | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0, 0.1]  | _    | "Spec maturity maint costs "
    κmat::F            | 0.05            | _               | Beta(2.0, 2.0)  | [0.0, 1.0]  | _    | "Reserve flux allocated to development/reprod."
    threshold::Mo      | 1.0             | mol             | Beta(2.0, 2.0)  | [1e-3, 20.0] | true | "Structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  
    # w_M::GMo         | 25.0            | g*mol^-1        | Beta(2.0, 2.0)  | [0.0, 1.0]  | _   | "Mol-weight of shoot maturity reserve:"
end

@columns struct Maturity{MoMoD,F,Mo} <: AbstractMaturity
    # Field            | Default         | Unit            | Prior           | Limits      | Log  | Description
    # n_N_M::MoMo      | 0.05            | mol*mol^-1      | Gamma(2.0, 2.0) | [0.0, 1.0]  | _    | "N/C use for maturity"
    j_E_mat_mai::MoMoD | 0.001           | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0, 0.1]  | _    | "Spec maturity maint costs "
    κmat::F            | 0.05            | _               | Beta(2.0, 2.0)  | [0.0, 1.0]  | _    | "Reserve flux allocated to development/reprod."
    threshold::Mo      | 1.0             | mol             | Beta(2.0, 2.0)  | [1e-3, 20.0] | true | "Structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  
    # w_M::GMo         | 25.0            | g*mol^-1        | Beta(2.0, 2.0)  | [0.0, 1.0]  | _   | "Mol-weight of shoot maturity reserve:"
end


"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Nothing, o, u) = nothing

maturity!(f::Maturity, o, u) = begin
    seedset = u.M > f.threshold ? u.M : zero(u.M) 
    mat = κmat(f) * flux1(o)[:E,:ctb]
    flux(o)[:M,:gro] = mat - seedset * unit(mat)/unit(seedset)
    mat_mai = f.j_E_mat_mai * tempcorrection(o) * u.V
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
    # reserve_loss!(o, mat_mai)
    # conversion_loss!(o, mat, f.n_N_M)
end

κmat(maturity_pars::Maturity) = maturity_pars.κmat
