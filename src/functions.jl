
"""
    half_saturation(max, half, x)

Half satration curve.
"""
half_saturation(max, half, x) = max/(oneunit(half/x) + half/x)
