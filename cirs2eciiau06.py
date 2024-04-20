from iau06pna import iau06pna
from iau06pnb import iau06pnb
from iau06xys import iau06xys

def cirs2eciiau06(rcirs, vcirs, acirs, ttt, option, ddx, ddy):
    """
    This function transforms a vector from the CIRS frame to
    the ECI mean equator mean equinox (GCRF).

    Author: David Vallado, 719-573-2600, 2 May 2020

    Inputs:
        rcirs: position vector earth fixed (km)
        vcirs: velocity vector earth fixed (km/s)
        acirs: acceleration vector earth fixed (km/s^2)
        ttt: Julian centuries of TT (centuries)
        option: which approach to use ('a' - 2000a, 'b' - 2000b, 'c' - 2000xys)
        ddx: EOP correction for x (rad)
        ddy: EOP correction for y (rad)

    Outputs:
        reci: position vector ECI (km)
        veci: velocity vector ECI (km/s)
        aeci: acceleration vector ECI (km/s^2)

    Locals:
        nut: transformation matrix for ire-gcrf

    Coupling:
        iau00era: rotation for Earth rotation pef-ire
        iau00xys: rotation for prec/nut ire-gcrf

    References:
        Vallado, 2004, 205-219
    """

    # ---- ceo based, iau2006
    if option == 'c':
        x, y, s, pnb = iau06xys(ttt, ddx, ddy)

    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
        lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = iau06pna(ttt)

    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
        lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = iau06pnb(ttt)

    # ---- perform transformations
    reci = pnb @ rcirs
    veci = pnb @ vcirs
    aeci = pnb @ acirs

    return reci, veci, aeci
