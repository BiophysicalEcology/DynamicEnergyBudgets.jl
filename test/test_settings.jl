using DynamicEnergyBudgets

# parameters for plant
# The values of these parameters are very provisional

function test_settings(tspan::Tspan; save_intermediate=true, 
                       use_environment = false, state_type = StatePVCNE, 
                       environment = DataFrame(), years=1, 
                       timestep_days=1.0/(24.0*60*60), location=[148, -32]) # Everything was per day despite units in seconds and days.
    @deb_settings quote
        settings = begin
            u0 = u0
            tspan = tspan
            environment = environment
            use_environment = use_environment
            apply_environment! = f(x...) = nothing
            save_intermediate = save_intermediate
            timestep_days = timestep_days
            state_type = state_type
        end
        structures = begin
            shoot = begin 
                functions = begin
                    assim = shoot_assimilation!
                    assim_sub = photosynthesis_kooijman
                    area = area_mass_kooijman
                    rate = find_rate
                end
                params = begin
                    J_L_F = 5.0, "mol.mol⁻¹s⁻¹, flux of useful photons"
                    J_L_K = 1.0, "mol.mol⁻¹s⁻¹, half-saturation flux of useful photons"
                    j_L_Amax = 50.0, "mol.mol⁻¹.sec⁻¹, max spec uptake of useful photons"
                    j_C_Amax = 50.0, "mol.mol⁻¹.sec⁻¹, max spec uptake of carbon dioxide"
                    j_O_Amax = 0.0001, "mol.mol⁻¹.sec⁻¹, max spec uptake of oxygen"
                    k_C_binding = 1.0, "mols.s⁻¹, scaling rate for carbon dioxide"
                    k_O_binding = 1.0, "mols.s⁻¹, scaling rate for oxygen"
                    k_E = 0.2, "mol.d⁻¹, shoots reserve turnover rate"
                    k_EC = 0.2, "mol.d⁻¹, shoots C-reserve turnover rate"
                    k_EN = 0.2, "mol.d⁻¹, shoots N-reserve turnover rate"
                    j_E_mai = 0.001, "mol.mol⁻¹.d⁻¹, shoots spec somatic maint costs."
                    j_E_rep_mai = 0.001, "mol.mol⁻¹.d⁻¹, shoots spec maturity maint costs"
                    j_P_mai = 0.01, "mol.mol⁻¹.d⁻¹ shoot product formation linked to maintenance"
                    K_C = 1.0, "mol/L, half-saturation concentration of carbon dioxide"
                    K_O = 1.0, "mol/L, half-saturation concentration of oxygen"
                    X_C = 10.0, "mol/L carbon dioxide"
                    X_O = 100.0, "mol/L oxygen"
                    M_Vgerm = 0.5, "mol, shoots structural mass at germination"
                    M_Vrep = 10.0, "mol, shoots structural mass at start reproduction"
                    M_Vref = 1.0, "mol"
                    M_Vscaling = 100.0, "mol shoots scaling mass"
                    y_E_CH_NO = 1.5, "mol/mol, from shoots C-reserve to reserve, using nitrate"
                    y_V_E = 0.7, "mol/mol, from shoots reserve to structure"
                    y_E_ET = 0.8, "mol/mol, from shoots reserve to roots reserve"
                    y_EN_ENT = 1.0, "mol/mol, from roots N-reserve to shoots N-reserve"
                    y_E_EN = 0.5, "mol/mol, from shoots N-reserve to reserve"
                    y_P_V = 0.02, "mol/mol, shoot product formation linked to growth"
                    κEC = 0.2, "-, shoots non-processed C-reserve returned to C-reserve, the remaining fraction is translocated to the root"
                    κEN = 0.5, "-, shoots non-processed N-reserve returned to N-reserve"
                    κsoma = 0.6, "-, shoots reserve flux allocated to soma"
                    κrep = 0.05, "-, shoots reserve flux allocated to development/reprod."
                    n_N_P = 0.0, "- N/C in shoot product (wood)"
                    n_N_V = 0.15, "- N/C in shoot structure"
                    n_N_EC = 0.0, "- N/C in shoot C-reserve"
                    n_N_EN = 10.0, "- N/C in shoot N-reserve"
                    n_N_E = 0.2, "- N/C in shoot reserve"
                    w_P = 25.0, "g/mol, mol-weight of shoot product (wood)"
                    w_V = 25.0, "g/mol, mol-weight of shoot structure"
                    w_EC = 25.0, "g/mol, mol-weight of shoot C-reserve"
                    w_EN = 25.0, "g/mol, mol-weight of shoot N-reserve"
                    w_E = 25.0, "g/mol, mol-weight of shoot reserve"

                    REFERENCE_TEMP = 310.0, (273.0, 325.0), "K, temp for which rate pars are given"
                    ARRH_TEMP = 2000.0, (200.0, 4000.0), "K, Arrhenius temp"
                    LOWER_BOUNDARY = 280.0, (273.0, 325.0), "K, lower boundary tolerance range"
                    ARRH_LOWER = 20000.0, (2000.0, 40000.0), "K, Arrhenius temp for lower boundary"
                    UPPER_BOUNDARY = 315.0, (273.0, 325.0), "K, upper boundary tolerance range"
                    ARRH_UPPER = 70000.0, (7000.0, 140000.0), "K Arrhenius temp for upper boundary"
                    K_autophagy = 0.000001, (0.0000001, 0.00001)
                end
            end
            root = begin
                functions = begin
                    assim = root_assimilation!
                    assim_sub = nitrogen_uptake
                    area = area_mass_kooijman
                    rate = find_rate
                end
                params = begin
                    j_NH_Amax = 0.5, "mol.mol⁻¹.s⁻¹, max spec uptake of ammonia"
                    j_NO_Amax = 0.5, "mol.mol⁻¹.s⁻¹, max spec uptake of nitrate"
                    j_E_mai = 0.003, "mol.mol⁻¹.d⁻¹, roots spec somatic maint costs"
                    j_E_rep_mai = 0.0, "mol.mol⁻¹.d⁻¹, roots spec maturity maint costs"
                    j_P_mai = 0.01, "mol.mol⁻¹.d⁻¹, root product formation linked to maintenance"
                    K_NH = 10.0, "mol/L, half-saturation concentration of ammonia"
                    K_NO = 10.0, "mol/L, half-saturation concentration of nitrate"
                    K_H = 1.0, "mol/L, half-saturation concentration of water"
                    X_NH = 5.0, "mol/L ammonia"
                    X_NO = 10.0, "mol/L concentration of nitrate see e.g. [@crawford1998molecular]"
                    X_H = 10.0, "mol/L"
                    n_N_P = 0.0, "- N/C in root product (wood)"
                    n_N_V = 0.15, "- N/C in root structure"
                    n_N_EC = 0.0, "- N/C in root C-reserve"
                    n_N_EN = 10.0, "-, N/C in root N-reserve"
                    n_N_E = 0.2, "-, N/C introot reserve"
                    y_E_CH_NO = 1.5, "mol/mol, from roots C-reserve to reserve, using nitrate"
                    y_E_CH_NH = 1.25, "mol/mol, from roots C-reserve to reserve, using ammonia"
                    y_EC_ECT = 1.0, "mol/mol, from shoots C-reserve to roots C-reserve"
                    y_E_EN = 0.3, "mol/mol, from roots N-reserve to reserve"
                    κEC = 0.5, "-, roots  non-processed C-reserve returned to C-reserve, the remaining fraction is translocated to the shoot"
                    κEN = 0.2, "-, roots  non-processed N-reserve returned to N-reserve"
                    κsoma = 0.5, "-, roots  reserve flux allocated to soma"
                    κrep = 0.0, "-, shoots reserve flux allocated to development/reprod."
                    M_Vgerm = 0.3, "mol, roots structural mass at germination"
                    M_Vrep = 10.0, "mol, roots structural mass at start reproduction"
                    ρNO = 0.7, "-, weights preference for nitrate relative to ammonia."
                end
            end
        end
    end
end
