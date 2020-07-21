const parconv = 4.57mol*W^-1*s^-1

"""
    ManualTemperature(airtemperature, soiltemperature)

Environment for simple manual temperature control.
"""
struct ManualTemperature{T}
    airtemperature::T
    soiltemperature::T
end

"""
    apply_environment!(organism::Organism, organs::Tuple, u, t) 

Apply environmental conditions to an organism. 
"""
apply_environment!(organism::AbstractOrganism, organs, u, t) = 
    apply_environment!(environment(organism), organism, organs, u, t)

"""
    apply_environment!(organism, env::Nothing, organs, u, t)

No environment included, do nothing. Initial environment variables
will be used for the whole simulation
"""
apply_environment!(env::Nothing, organism, organs, u, t) = nothing

"""
    apply_environment!(plant, env::ManualTemperature, organs, u, t)

Apply an environment to a plant using live user input for variables,
instead of data.
"""
apply_environment!(env::ManualTemperature, plant::Plant, organs, u, t) = begin
    shoot, root = organs
    update_temp!(shoot, env.airtemperature)
    update_temp!(root, env.soiltemperature)
end

"""
    update_temp!(o, temp)

Set the temperature, and update the temperature correction variable.
"""
update_temp!(o, temp) = begin
    set_temp!(o, temp)
    set_tempcorrection!(o, tempcorr(tempcorr_pars(o), temp))
end
