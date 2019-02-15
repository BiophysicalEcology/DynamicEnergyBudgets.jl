abstract type AbstractMaturity end

" Maturity parameters. Seperated to make maturity modeling optional, reducing complexity "
@columns struct Maturity{MoMo,MoMoD,F,Mo} <: AbstractMaturity
    # Field            | Default         | Unit            | Prior           | Limits      | Log | Description
    n_N_M::MoMo        | 0.1             | mol*mol^-1      | Gamma(2.0, 2.0) | [0.0, 1.0]  | _   | "N/C use for maturity"
    j_E_mat_mai::MoMoD | 0.001           | mol*mol^-1*d^-1 | Beta(2.0, 2.0)  | [0.0, 0.1]  | _   | "Shoots spec maturity maint costs "
    κmat::F            | 0.05            | _               | Beta(2.0, 2.0)  | [0.0, 1.0]  | _   | "Shoots reserve flux allocated to development/reprod."
    M_Vmat::Mo         | 10.0            | mol             | Beta(2.0, 2.0)  | [0.0, 20.0] | _   | "Shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?  w_M::GMo           | 25.0            | g*mol^-1        | Beta(2.0, 2.0)  | [0.0, 1.0]   | "Mol-weight of shoot maturity reserve:"
end


"""
    maturity!(f, o, u)
Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Nothing, o, u) = nothing

maturity!(f::Maturity, o, u) = begin
    flux(o)[:M,:gro] = mat = f.κmat * flux1(o)[:E,:ctb]
    mat_mai = f.j_E_mat_mai * tempcorrection(v) * u.V # min(u[:V], f.M_Vmat))
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
    # reserve_loss!(o, mat_mai)
    # conversion_loss!(o, mat, f.n_N_M)
end
