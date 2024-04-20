from cubicspl import cubicspl
from cubic import cubic

def cubicinterp(p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein):
    """
    This function performs a cubic spline. Four points are needed.

    Author: David Vallado, 719-573-2600, 1 Dec 2005

    Parameters:
        p1a, p1b, p1c, p1d (float): Points for cubic spline 1.
        p2a, p2b, p2c, p2d (float): Points for cubic spline 2.
        valuein (float): Input value.

    Returns:
        float: Interpolated value.

    References:
        Vallado 2013, 1027
    """

    # -------- assign function points ---------
    ac0, ac1, ac2, ac3 = cubicspl(p1a, p1b, p1c, p1d)
    kc0, kc1, kc2, kc3 = cubicspl(p2a, p2b, p2c, p2d)

    # recover the original function values
    # use the normalized time first, but at an arbitrary interval
    r1r, r1i, r2r, r2i, r3r, r3i = cubic(kc3, kc2, kc1, kc0 - valuein, 'R')
    # print('cubic', ac0, ac1, kc0, r1r, r2r)
    if -0.000001 <= r1r <= 1.001:
        value = r1r
    elif -0.000001 <= r2r <= 1.001:
        value = r2r
    elif -0.000001 <= r3r <= 1.001:
        value = r3r
    else:
        value = 0.0
        print(f"error in cubicinterp root {valuein} {r1r} {r2r} {r3r}")

    answer = ac3 * value ** 3 + ac2 * value ** 2 + ac1 * value + ac0
    return answer
