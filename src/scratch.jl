using Revise
using MechanisticModels
using DynamicEnergyBudgets
using AutoInteract
using JLD
using Reactive
using Interact
using ProfileView
# gr()

include("../DynamicEnergyBudgets/src/settings.jl")

use_environment = true
save_intermediate = true
environment = []
# environment = DynamicEnergyBudgets.load_nichemap("Northcote, victoria", 1)
# save("environment.jld", "environment", environment)
environment = load("environment.jld")["environment"]
settings = build_settings((0.0, 2000.0), save_intermediate=save_intermediate, use_environment=use_environment, environment=environment)
MechanisticModels.cascade_update_params!(settings.structures)
apply(DynamicEnergyBudgets.scale_time_dependent_params!, settings.structures, settings.timestep_days)
apply(DynamicEnergyBudgets.initialise_params!, settings.structures)
# plot(environment[:soil])

w = deb_widgets(settings)
p = deb_plottables(settings)
sp = get_signals(p); sp.value
sw = get_signals(w); sw.value
i = make_interface(p)
delete!(p, :structures)
settings = apply_all(settings, sw.value)
MechanisticModels.cascade_update_params!(settings.structures)
apply(DynamicEnergyBudgets.scale_time_dependent_params!, settings.structures, settings.timestep_days)
apply(DynamicEnergyBudgets.initialise_params!, settings.structures)

plot_all(settings, sp.value, 1:10) 

struct x
    n::typeof(u"molN")
end

@time for i = 1:100 integrate(settings) end; nothing

prof(view = true)

function prof(; maxdepth = 40, view = false)
    Profile.clear()
    integrate(settings)
    for i = 1:10
        @profile m = integrate(settings)
    end
    !view || ProfileView.view();
    Profile.print(maxdepth = maxdepth)
end

# function shape(settings; tmax = 100)
#   pyplot()
#   (states, fluxes) = integrate_custom(settings, tmax)
#   verts = [(-1.0,1.0),(-1.0,-1.0),(1.0,-1.0),(1.0,1.0)]
#   sizesS = states[:, VS]
#   sizesR = states[:, VR]
#   println("Animating...")
#   labels = ["Shoot Structure" "Root Structure" "Reserve"]
#   anim = @animate for i = 1:100:tmax
#     x = [1.0]
#     y = [1.0]
#     an_x = [2.0]
#     marker = :circle, (states[i, VS] + states[i, ES])/10
#     plot(x, y, marker = marker, color = :green,
#          bg = :grey, fg = :black,
#          xlim = (0,3), ylim = (-2,2),
#          leg = false,
#          label = labels,
#          annotations = (an_x, y, "Shoot volume with\nreserve portion")
#         )
#     marker = :circle, (states[i, ES] / 10)
#     plot!(x, y, marker = marker, color = :purple)
#     y = [-1.0]
#     marker = :circle, (states[i, VR] + states[i, ER]) / 20
#     plot!(x, y, marker = marker, color = :brown,
#           annotations = (an_x, y, "Root volume with\nreserve portion\n(so much reserve!)")
#          )
#     marker = :circle, (states[i, ER]/20)
#     plot!(x, y, marker = marker, color = :purple)
#   end
#   return anim
# end
