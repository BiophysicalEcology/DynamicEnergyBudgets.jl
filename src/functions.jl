"""
    non_growth_flux(ureserve, turnover, r)
Returns the current non_growth_flux flux at rate r,
or the flux as a proportion of u[V], depending on ureserve values.
"""
non_growth_flux(reserve, turnover, r) = reserve * (turnover - r)

"""
    half_saturation(max, half, x)
Half satration curve.
"""
half_saturation(max, half, x) = max/(oneunit(half/x) + half/x)
