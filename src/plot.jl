import Plots.px

function debplot(settings, signals, plottable)
    # settings.u0 = collect(Iterators.flatten([s[:u] for s in signals[:structures]]))
    settings = apply_all(settings, signals)
    structures = settings.structures
    MechanisticModels.cascade_update_params!(structures)
    apply(DynamicEnergyBudgets.scale_time_dependent_params!, structures, settings.timestep_days)
    apply(DynamicEnergyBudgets.initialise_params!, structures)
    tmin, tmax = floor.(Int, settings.tspan) .+ 1
    sol = integrate(settings)
    stateplot = plot(sol, xlab="", link=:x, plotdensity=200, margin=0px, linewidth=3, alpha=0.7)
    selected = plot_all(settings, plottable, tmin:tmax)
    if selected != nothing
        l = Plots.GridLayout(2, 1)
        return plot(stateplot, selected, layout=l, size=(1300, 700))
    else
        return plot(stateplot, size=(1300, 700))
    end
end

function plot_function(func)
    temp_keys = [:REFERENCE_TEMP, :ARRH_TEMP, :LOWER_BOUNDARY, :ARRH_LOWER, :UPPER_BOUNDARY, :ARRH_UPPER]
    area_keys = [:M_Vref, :M_Vscaling]
    colors = distinguishable_colors(10)
    state_toggles = [[true, true, true, true, true, true],
                     [true, true, true, true, true, true]]  
    (state_vars, state_labels, state_colors) = select_state_vars(structures, state_toggles)
    temp_args = (temp_keys, tempcorr, 270.0:2.0:330.0)
    area_args = (area_keys, area_mass_kooijman, 0.0:2.0:400.0)
    num_subs = 0
    for args in (temp_args, area_args)
        num_subs += 1
        push_func_subplot(structures, args...)
    end
    if num_subs > 0
        height = num_subs < 2 ? 0.5pct : :auto
        subplots_layout = Plots.GridLayout(num_subs, 1, width=0.2pct)
        col_layout = Plots.GridLayout(1, 2)
        col_layout[1, 1] = row_layout
        col_layout[1, 2] = subplots_layout
        master_layout = col_layout
    end
end

function select_state_vars(structures, state_toggles)
    num_structures = length(state_toggles)
    len = length(state_toggles[1])
    scols = distinguishable_colors(len)
    state_labels = Array{String,2}(num_structures,len)
    state_colors = Array{RGB,2}(num_structures,len)
    state_vars = Array{Tuple,2}(num_structures,len)
    state_names = fieldnames(structures[1].u)
    for (n, st) in enumerate(state_toggles)
        i = 0
        for (index, sig) in enumerate(st)
            if sig
                i += 1
                state_vars[n, i] = (0, index)
                state_labels[n, i] = string(state_names[index])
                state_colors[n, i] = scols[index]
            end
        end
    end
    return (state_vars, state_labels, state_colors)
end

function func_subplot(structures, sub_keys, func, range)
    p = plot(legend=:none, linewidth=3, alpha=0.7)
    for s in structures
        sub_args = [getfield(s.vars, key) for key in sub_keys] 
        plot!(p, x -> func(x, sub_args...), range)
    end
    return p
end
