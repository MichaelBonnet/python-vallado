def recovpar(p1, p2, p3, p4, root):
    """
    Recover the time and function values in parabolic blending routines.

    Args:
    p1, p2, p3, p4 (float): Function values used for blending.
    root (float): Root used as variable.

    Returns:
    float: Function value.
    """
    # Set up function from C-39
    # acut3*x**3 + acut2*x**2 + etc
    acut0 = p2
    acut1 = (-p1 + p3) * 0.5
    acut2 = p1 - 2.5 * p2 + 2.0 * p3 - 0.5 * p4
    acut3 = -0.5 * p1 + 1.5 * p2 - 1.5 * p3 + 0.5 * p4

    # Recover the variable value
    funvalue = acut3 * root**3 + acut2 * root**2 + acut1 * root + acut0

    return funvalue
