using AutoInteract

function deb_widgets(settings)
    widgets = make_widgets(settings)
    for s in widgets[:structures]
        delete!.(s, [:A, :params, :init_params, :param_ids, :J1, :J, :flags, :functions, :assim_state])
        # Put parameter sliders last
        param_specs2 = s[:param_specs]
        pop!(s, :param_specs)
        s[:param_specs] = param_specs2
    end
    delete!.(widgets, [:u0, :environment])
    widgets[:tspan] = (make_widgets(settings.tspan[1], "tmin", 0.0:1.0:settings.tspan[2]),
                       make_widgets(settings.tspan[2], "tmax", 0.0:1.0:settings.tspan[2]))
    return widgets
end

function deb_plottables(settings)
    widgets = make_plottables(settings)
    for s in widgets[:structures]
        delete!.(s, [:A, :params_specs, :params, :init_params, :param_ids, :J1, :J, :flags, :functions, :assim_state])
    end
    # if :environment in keys(widgets)
    #   for k in [:nyears, :timeinterval, :minshade, :maxshade, :dim, :ALTT, :REFL, :longlat, :DEP, :MAXSHADES]
    #     delete!(widgets[:environment], k)
    #   end
    # end
    delete!(widgets, :u0)
    return widgets
end

