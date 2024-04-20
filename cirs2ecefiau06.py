import numpy as np
from iau06era import iau06era
from iau06gst import iau06gst
from iau06pna import iau06pna
from iau06pnb import iau06pnb
from iau06xys import iau06xys
from fundarg import fundarg
from polarm import polarm
from constastro import earthrot

def cirs2ecefiau06(rcirs, vcirs, acirs, ttt, jdut1, lod, xp, yp, option, ddx, ddy):
    """
    This function transforms a vector from the CIRS (GCRF), to an Earth fixed (ITRF) frame.
    The results take into account the effects of sidereal time, and polar motion.

    Author: David Vallado, 719-573-2600, 2 May 2020

    Revisions:

    Inputs:
        rcirs: position vector CIRS (km)
        vcirs: velocity vector CIRS (km/s)
        acirs: acceleration vector CIRS (km/s^2)
        ttt: Julian centuries of TT (centuries)
        jdut1: Julian date of UT1 (days from 4713 BC)
        lod: Excess length of day (sec)
        xp: Polar motion coefficient (arc sec)
        yp: Polar motion coefficient (arc sec)
        option: Which approach to use (a-2000a, b-2000b, c-2000xys)

    Outputs:
        recef: position vector Earth fixed (km)
        vecef: velocity vector Earth fixed (km/s)
        aecef: acceleration vector Earth fixed (km/s^2)
    """

    # ---- ceo based, iau2000
    if option == 'c':
        x, y, s, pnb = iau06xys(ttt, ddx, ddy)
        st = iau06era(jdut1)

    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = iau06pna(ttt)
        gst, st = iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate)

    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = iau06pnb(ttt)
        gst, st = iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate)

    pm = polarm(xp, yp, ttt, '06')

    # Setup parameters for velocity transformations
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0, 0, thetasa])

    rpef = np.dot(st.T, rcirs)
    recef = np.dot(pm.T, rpef)

    vpef = np.dot(st.T, vcirs) - np.cross(omegaearth, rpef)
    vecef = np.dot(pm.T, vpef)

    temp = np.cross(omegaearth, rpef)
    aecef = np.dot(pm.T, (np.dot(st.T, acirs) - np.cross(omegaearth, temp) - 2.0 * np.cross(omegaearth, vpef)))

    return recef, vecef, aecef
