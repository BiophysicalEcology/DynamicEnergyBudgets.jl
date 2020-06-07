abstract type AbstractDEBCore end

# TODO generalise this.

@columns struct DEBCore{MoMo,GMo} <: AbstractDEBCore
    # Field      | Default | Unit     | Limits       | Log   | Description
    y_V_E::MoMo  | 0.7     | _        | [0.0, 1.0]   | _     | "Yield from reserve to structure"
    y_E_EC::MoMo | 0.7     | _        | [1e-6, 1.0]  | false | "Yield from C-reserve to general reserve"
    y_E_EN::MoMo | 30.0    | _        | [1.0, 50.0]  | false | "Yeild from N-reserve to general reserve"
    n_N_V::MoMo  | 0.03    | _        | [0.0, 0.1]   | _     | "Nitrogen per Carbon in structure"
    n_N_E::MoMo  | 0.025   | _        | [0.0, 0.1]   | _     | "Nitrogen per Carbon in reserve"
    w_V::GMo     | 25.0    | g*mol^-1 | [15.0, 40.0] | _     | "Mol-weight of shoot structure"
    # w_N::GMo   | 25.0    | g*mol^-1 | [15.0, 40.0] | _     | "Mol-weight of shoot N-reserve"
    # w_C::GMo   | 25.0    | g*mol^-1 | [12.0, 40.0] | _     | "Mol-weight of shoot C-reserve"
    # w_E::GMo   | 25.0    | g*mol^-1 | [15.0, 40.0] | _     | "Mol-weight of shoot reserve"
end

