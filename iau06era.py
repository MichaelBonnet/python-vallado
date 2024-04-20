from math import cos, sin, pi

def iau06era(jdut1):
    """
    Calculate the transformation matrix that accounts for the effects of sidereal time via the Earth rotation angle.

    Parameters:
    - jdut1: Julian date of UT1 (days)

    Returns:
    - st: Transformation matrix for PEF-IRE
    """



    # Julian centuries of UT1
    tut1d = jdut1 - 2451545.0

    # Earth rotation angle
    era = (2 * pi) * (0.7790572732640 + 1.00273781191135448 * tut1d)
    era = era % (2 * pi)

    # Print ERA if needed
    # if iauhelp == 'y':
    #     print(f'era{era * 180 / pi:11.7f}')

    # Transformation matrix
    st = [[cos(era), -sin(era), 0.0],
          [sin(era),  cos(era), 0.0],
          [0.0,       0.0,      1.0]]

    return st
