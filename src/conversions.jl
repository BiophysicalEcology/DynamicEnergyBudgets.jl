export watts_to_light_mol, light_mol_to_watts, water_content_to_mols_per_litre,
       celcius_to_kelvin, kelvin_to_celcius, second_to_day, second_to_hour,
       hour_to_day, mmol_to_mol, μmol_to_mol, gram_to_mol,
       fraction_per_litre_gas_to_mols

# Return mol/m²s. Usually this is measured in μmol.
watts_to_light_mol(watts::Number)::Float64 = watts * 4.57e-6
light_mol_to_watts(light_mol::Number)::Float64 = light_mol / 4.57e-6

water_content_to_mols_per_litre(wc::Number)::Float64 = wc * 55.5 # L/L of water to mol/L 
celcius_to_kelvin(degrees_C::Number)::Float64 = degrees_C + 273.15
kelvin_to_celcius(K::Number)::Float64 = K - 273.15

second_to_day(second::Number)::Number = second * 60 * 60 * 24
second_to_hour(second::Number)::Number = second * 60 * 60
hour_to_day(hour::Number)::Number = hour * 24

mmol_to_mol(mmol::Number) = mmol * 1e-3
μmol_to_mol(μmol::Number) = μmol * 1e-6

gram_to_mol(mol_weight::Number, g::Number) = g / mol_weight
fraction_per_litre_gas_to_mols(frac::Number)::Number = frac / 22.4
