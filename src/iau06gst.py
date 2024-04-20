import numpy as np
from iau06in import iau06in

def iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega,
             lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate):
    """
    This function finds the IAU2006 Greenwich Sidereal Time.

    Author: David Vallado, 719-573-2600, 16 Jul 2004

    Inputs:
        jdut1: Julian date of UT1 (days from 4713 BC)
        ttt: Julian centuries of TT
        deltapsi: Change in longitude (rad)
        l, l1, f, d, omega: Delaunay elements (rad)
        lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep: Planetary values (rad)
        precrate: Precession rate (rad)

    Outputs:
        gst: Greenwich Sidereal Time (0 to twopi rad)
        st: Transformation matrix

    Locals:
        temp: Temporary variable for reals (rad)
        tut1d: Days from the Jan 1, 2000 12 h epoch (UT1)

    References:
        Vallado, 2004, 216
    """
    # Constants
    deg2rad = np.pi / 180.0
    convrt = np.pi / (180.0 * 3600.0)

    # Initialize arrays
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()

    # Precompute powers of ttt
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2

    # Mean obliquity of the ecliptic
    epsa = 84381.406 - 46.836769 * ttt - 0.0001831 * ttt2 + 0.00200340 * ttt3 - 0.000000576 * ttt4 \
           - 0.0000000434 * ttt5  # arcseconds
    epsa = (epsa / 3600.0) % 360.0  # degrees
    epsa = epsa * deg2rad  # radians

    # Evaluate the ee complementary terms
    gstsum0 = 0.0
    for i in range(33, 0, -1):
        tempval = agsti[i-1, 0] * l + agsti[i-1, 1] * l1 + agsti[i-1, 2] * f + agsti[i-1, 3] * d + agsti[i-1, 4] * omega + \
                  agsti[i-1, 5] * lonmer + agsti[i-1, 6] * lonven + agsti[i-1, 7] * lonear + agsti[i-1, 8] * lonmar + \
                  agsti[i-1, 9] * lonjup + agsti[i-1, 10] * lonsat + agsti[i-1, 11] * lonurn + agsti[i-1, 12] * lonnep + \
                  agsti[i-1, 13] * precrate
        gstsum0 += agst[i-1, 0] * np.sin(tempval) + agst[i-1, 1] * np.cos(tempval)  # rad
    gstsum1 = 0.0
    for j in range(1, -1, -1):
        i = 33 + j
        tempval = agsti[i-1, 0] * l + agsti[i-1, 1] * l1 + agsti[i-1, 2] * f + agsti[i-1, 3] * d + agsti[i-1, 4] * omega + \
                  agsti[i-1, 5] * lonmer + agsti[i-1, 6] * lonven + agsti[i-1, 7] * lonear + agsti[i-1, 8] * lonmar + \
                  agsti[i-1, 9] * lonjup + agsti[i-1, 10] * lonsat + agsti[i-1, 11] * lonurn + agsti[i-1, 12] * lonnep + \
                  agsti[i-1, 13] * precrate
        gstsum1 += agst[i-1, 0] * ttt * np.sin(tempval) + agst[i-1, 1] * ttt * np.cos(tempval)

    eect2000 = gstsum0 + gstsum1 * ttt  # rad

    # Equation of the equinoxes
    ee2000 = deltapsi * np.cos(epsa) + eect2000  # rad

    # Earth rotation angle
    tut1d = jdut1 - 2451545.0
    era = 2.0 * np.pi * (0.7790572732640 + 1.00273781191135448 * tut1d)
    era = era % (2.0 * np.pi)  # rad

    # Greenwich mean sidereal time, IAU 2000
    gmst2000 = era + (0.014506 + 4612.156534 * ttt + 1.3915817 * ttt2 \
                      - 0.00000044 * ttt3 + 0.000029956 * ttt4 + 0.0000000368 * ttt5) * convrt  # " to rad

    gst = gmst2000 + ee2000  # rad

    # Transformation matrix
    st = np.zeros((3, 3))
    st[0, 0] = np.cos(gst)
    st[0, 1] = -np.sin(gst)
    st[0, 2] = 0.0

    st[1, 0] = np.sin(gst)
    st[1, 1] = np.cos(gst)
    st[1, 2] = 0.0

    st[2, 0] = 0.0
    st[2, 1] = 0.0
    st[2, 2] = 1.0

    return gst, st
