export watts_to_light_mol, light_mol_to_watts, water_content_to_mols_per_litre,
       fraction_per_litre_gas_to_mols

watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L 
fraction_per_litre_gas_to_mols(frac) = frac / 22.4
