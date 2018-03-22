using DynamicEnergyBudgets
using NicheMap
using Unitful

function build_settings(tspan; environment=[], use_environment=false, 
                        save_intermediate=false, timestep_days=1.0/24.0) 
    state_type = StatePVMCNE
    if state_type == StatePVMCN
        u0 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 10.0] # Initial value
        # u0 = [0.0, 1e-4, 1e-2, 1e-2, 0.0, 1e-4, 8.0, 2.0, 0.0, 1e-4, 8.0, 2.0] # Initial value
    elseif state_type == StatePVCNE
        u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0] # Initial value
        # u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0, 0.0, 1e-4, 1e-4, 1e-4, 10.0] # Initial value
    elseif state_type == StatePVMCNE
        u0 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0] # Initial value
        # u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0, 0.0, 1e-4, 1e-4, 1e-4, 10.0] # Initial value
    end

    @deb_settings quote
        settings = begin
            u0 = u0
            tspan = tspan
            environment = environment
            use_environment = use_environment
            apply_environment! = apply_energy_balance!
            save_intermediate = save_intermediate
            timestep_days = timestep_days
            state_type = state_type
        end
        structures = begin
            Leaf = begin
                functions = begin
                    area = area_mass_kooijman
                    assim = shoot_assimilation!
                    assim_sub = photomoduleC3
                    rate = find_rate
                end
                params = begin 
                    # J_L_F = watts_to_light_mol(800.0) # mol/m²s, flux of useful photons
                    # TODO: should this acclimatise over time? this would require another state variable. 

                    # J_L_F should feed back to affect plant state: photosynthesis as a sensory 
                    # as well as energetic process [_@huner1998energy
                    # J_L_K = watts_to_light_mol(300.0), "mol/m²s, half-saturation flux of useful photons"
                    # Max specific uptake parameters that relate uptake to active surface area
                    # These are usually measured in μmol m⁻² s⁻¹ (wikipedia:Photosynthetic capacity).
                    # From [_@walker2014relationship :...full range of photosynthetically active radiation 
                    # (PAR 0–1500 μmol·m−2·s−1) and three levels of Vcmax (25, 50 & 90 μmol·m−2·s−1)
                    # j_L_Amax = μmol_to_mol(20.0) #[@2006seasonality] umol.m⁻².s⁻¹, max spec uptake of useful photons
                    # j_C_Amax = μmol_to_mol(90.0), "mol.m⁻².s⁻¹, max spec uptake of carbon dioxide"
                    # j_O_Amax = μmol_to_mol(0.001), "mol.m⁻².s⁻¹, max spec uptake of oxygen"
                    # Binding rates of gases to quantify photo-respiration 
                    # Is this braodly similar accross plants?
                    # k_C_binding = 1.0, :temp, "mols.s⁻¹, scaling rate for carbon dioxide "
                    # k_O_binding = 1.0, :temp, "mols.s⁻¹, scaling rate for oxygen"
                    # Turnover. This controls metabolic rate along with area/mass relation and 
                    # current reserves. 
                    k_E = "molE/dy", 0.2, (0.0, 1.0), (:exposed, :time, :temp), "shoots reserve turnover rate"
                    k_EC = "molEC/dy", 0.2, (0.0, 1.0), (:exposed, :time, :temp), "shoots C-reserve turnover rate"
                    k_EN = "molEN/dy", 0.2, (0.0, 1.0), (:exposed, :time, :temp), "shoots N-reserve turnover rate"
                    # Specific maintenancs costs in terms of reserve fluxes. 
                    # Must be per day not s, or these numbers are crazy.
                    # - these costs are paid to maintain structural mass and reproductive maturity
                    j_E_mai = "molE/molV*dy^⁻¹", 0.001, (0.0, 0.01), (:exposed, :time, :temp), "shoots spec somatic maint costs."
                    j_E_rep_mai = "molE/molV*dy^⁻¹", 0.001, (0.0, 0.01), (:exposed, :time, :temp), "shoots spec maturity maint costs "
                    # Production parameters
                    # - these parameters play no dynamic role but can dominate weights
                    j_P_mai = "molP/molE*dy^⁻¹", 0.01, (0.0, 0.1), (:exposed, :time, :temp), "shoot product formation linked to maintenance"

                    SLA = "m^²/g", 9.10, (5.0, 30.0), "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"

                    K_autophagy = "", 0.000001, (0.0000001, 0.00001)
                    K_C = "molC/L", fraction_per_litre_gas_to_mols(40.0/1e6), "half-saturation concentration of carbon dioxide"
                    K_O = "molO/L", fraction_per_litre_gas_to_mols(0.0021), "half-saturation concentration of oxygen"
                    X_C = "molC/L", fraction_per_litre_gas_to_mols(400.0/1e6), "carbon dioxide @ 400ppm"
                    X_O = "molO/L", fraction_per_litre_gas_to_mols(0.21), "oxygen (21% volume in air) "
                    # Life stage parameters
                    M_Vgerm = "molV", 0.5, "shoots structural mass at germination"
                    M_Vrep = "molV", 10.0, "shoots structural mass at start reproduction" # TODO: isn't this variable/seasonally triggered?

                    # Parameters that link active surface area to structural mass
                    # - they describe the development through V1- iso- and V0-morphs
                    # TODO: estimate these from functional traints SLA/LMA and/or 
                    # stem/branch wood specific gravity and final height? 
                    # The curve function probably also needs replacing.
                    M_Vref = "molV", 4.0, (0.4, 20.0), :exposed
                    M_Vscaling = "molV", 400.0, (40.0, 2000.0), :exposed, "shoots scaling mass"

                    # Partitioning parameters: dimensionless fractions
                    κEC = "", 0.2, (0.0, 1.0), :exposed, "shoots non-processed C-reserve returned to C-reserve,"
                    #    the remaining fraction is translocated to the root
                    κEN = "", 0.5, (0.0, 1.0), (:exposed, :time), "shoots non-processed N-reserve returned to N-reserve"
                    κsoma = "", 0.6, (0.0, 1.0), :exposed, "shoots reserve flux allocated to soma              "
                    κrep = "", 0.05, (0.0, 1.0), :exposed, "shoots reserve flux allocated to development/reprod."

                    y_V_E = "molE/molV", 0.7, (:exposed, :time), "from shoots reserve to structure: 0.3 lost as CO2?"
                    y_P_V = "molE/molP", 0.02, :exposed, "shoot product formation linked to growth: where does this carbon come from?"
                    y_E_CH_NO = "molEC/molE", 1.5, :exposed, "from shoots C-reserve to reserve, using nitrate: 0.75 EC per E."
                    y_E_EN = "molEN/molE", 0.5, :exposed, "from shoots N-reserve to reserve: 2 EN per E"
                    y_E_ET = "molE/molE", 0.8, :exposed, "from shoots reserve to roots reserve:"
                    y_EN_ENT = "molEN/molEN", 1.0, :exposed, "from roots N-reserve to shoots N-reserve"

                    # Chemical indices: elemental_nitrogen/elemental_carbon
                    # - only n_N_EN and n_N_E play a dynamic role
                    # - the others are only used to evaluate the nitrogen balance
                    n_N_P = "molN/molC", 0.0, "N/C in product (wood)"
                    n_N_V = "molN/molC", 0.15, "N/C in structure. Shouldnt this be identical to the reserve?"
                    n_N_EC = "molN/molC", 0.0, "N/C in C-reserve"
                    n_N_EN = "molN/molC", 10.0, "N/C in N-reserve"
                    n_N_E = "molN/molC", 0.2, "N/C in reserve. This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5 "

                    # Parameters that link moles to grams (wet weight)
                    # These should be calculated from first principles, not preset.
                    w_P = "g/mol" ,25.0, "mol-weight of shoot product (wood): 1mol * 12g/mol = 12g"
                    w_V = "g/mol" ,25.0, "mol-weight of shoot structure: 0.15mol * 14g/mol + 1mol * 12g/mol = 14.1g"
                    w_EC = "g/mol" ,25.0, "mol-weight of shoot C-reserve: 1mol * 12g/mol = 12g"
                    w_EN = "g/mol" ,25.0, "mol-weight of shoot N-reserve: 10mol * 14g/mol + 1mol * 12g/mol = 152g"
                    w_E = "g/mol" ,25.0, "mol-weight of shoot reserve: 0.2mol * 14g/mol + 1mol * 12g/mol = 14.8g"

                    REFERENCE_TEMP = "K", 310.0, (273.0, 325.0), "temp for which rate pars are given"
                    ARRH_TEMP = "K", 2000.0, (200.0, 4000.0), "Arrhenius temp"
                    LOWER_BOUNDARY = "K", 280.0, (273.0, 325.0), "lower boundary tolerance range temp_params[_UPPER_BOUNDARY] = 318.0 # K, upper boundary tolerance range"
                    ARRH_LOWER = "K", 20000.0, (2000.0, 40000.0), "Arrhenius temp for lower boundary"
                    UPPER_BOUNDARY = "K", 315.0, (273.0, 325.0), "upper boundary tolerance range"
                    ARRH_UPPER = "K", 70000.0, (7000.0, 140000.0), "Arrhenius temp for upper boundary"

                    LEAF_D = "", 0.01, (0.0001, 0.1)
                    LEAF_W = "", 0.01, (0.0001, 0.1)
                    LEAF_EMISSIVITY = "", 0.5, (0.0, 1.0)

                    Vcm25 = "", 110.0
                    Jm25 = "", 184.0
                    Vpm25 = "", -1.0
                    TPU25 = "", 15.0
                    Rd25 = "", 1.0
                    θ = "", 0.7
                    EaVc = "", 64800.0
                    Eaj = "", 37000.0
                    Hj = "", 220000.0
                    Sj = "", 710.0
                    Hv = "", 219400.0
                    EaVp = "", -1.0
                    Sv = "", -1.0
                    Eap = "", 47100.0
                    Ear = "", 66400.0
                    g0 = "", 0.02
                    g1 = "", 10.0
                    stoma_ratio = "", 0.5
                    leaf_width = "m", 0.05
                    leaf_angfact = "", 1.0
                    photo_flux_density = "μmol*m^-2*s^-1)", 1700.0
                    Tair = "C", 25.0
                    CO2 = "", 370.0
                    rel_humidity = "", 65.0
                    wind = "m/s", 0.8
                    pressure = "", 1.0
                end
            end
            #= Stem = begin
                functions = begin
                    area = area_mass_kooijman
                    assim = f(_...) = nothing
                    assim_sub = f(_...) = nothing
                    rate = find_rate
                end
                params = begin
                    M_Vref = 4.0, (0.4, 20.0), :exposed, "mol"
                end
            end # =#
            Root = begin
                functions = begin
                    area = area_mass_kooijman
                    assim = root_assimilation!
                    assim_sub = nitrogen_uptake
                    rate = find_rate
                end
                params = begin
                    j_NH_Amax = "mol*m^-2*s^-1", μmol_to_mol(50.0,), μmol_to_mol.((0.1, 1000.0)), "max spec uptake of ammonia"
                    j_NO_Amax = "mol*m^-2*s^-1", μmol_to_mol(50.0,), μmol_to_mol.((0.1, 1000.0)), "max spec uptake of nitrate"
                    j_E_rep_mai = "mol*mol^-1*dy^-1", 0.0, "roots spec maturity maint costs "

                    # Water affects root nutrient uptake via the saturation parameters
                    K_NH = "molNH/L", mmol_to_mol(10.0), "half-saturation concentration of ammonia"
                    K_NO = "molNOH/L", mmol_to_mol(10.0), "half-saturation concentration of nitrate"
                    K_H = "molH20/L", 10.0, "half-saturation concentration of water"
                    X_NH = "mol/L", mmol_to_mol(5.0), "ammonia"
                    X_NO = "mol/L", mmol_to_mol(10.0), "concentration of nitrate see e.g. [_@crawford1998molecular]"
                    X_H = "mol/L", 10.0 

                    # Life stage parameters
                    M_Vgerm = "mol", 0.3, "roots structural mass at germination"
                    M_Vrep = "mol", 10.0, "roots structural mass at start reproduction"

                    κEC = "", 0.5, (0.0, 1.0), :exposed, "roots  non-processed C-reserve returned to C-reserve"
                    #    the remaining fraction is translocated to the shoot
                    κEN = "", 0.2, (0.0, 1.0), :exposed, "roots  non-processed N-reserve returned to N-reserve"
                    κsoma = "", 0.5, (0.0, 1.0), :exposed, "roots  reserve flux allocated to soma               "
                    κrep = "", 0.0, (0.0, 1.0), "shoots reserve flux allocated to development/reprod."

                    ρNO = "", 0.7, (0.0, 1.0), "weights preference for nitrate relative to ammonia." # 1 or less but why?

                    # Yield coefficients (see also production parameters)
                    y_E_CH_NO = "molC/molE", 1.5, "from roots C-reserve to reserve, using nitrate"
                    y_E_CH_NH = "molC/molE", 1.25, "from roots C-reserve to reserve, using ammonia"
                    y_EC_ECT = "molC/molC", 1.0, "from shoots C-reserve to roots C-reserve"
                    y_E_EN = "molN/molE", 0.3, "from roots N-reserve to reserve: why is this different to shoots?"
                end
            end
        end
    end
end
