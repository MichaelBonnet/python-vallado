import numpy as np
from iau06in import iau06in
from fundarg import fundarg
from precess import precess
from rotations import rot1mat, rot2mat, rot3mat

def iau06pnb(ttt):
    """
    This function calculates the transformation matrix that accounts for the
    effects of precession-nutation in the IAU2000B theory.

    Author: David Vallado (719-573-2600), 16 Jul 2004

    Revisions:
    - Consolidate with IAU 2000, 14 Feb 2005

    Inputs:
        ttt: Julian centuries of TT

    Outputs:
        deltapsi: Change in longitude (rad)
        pnb: Transformation matrix for IRE-GCRF
        nut: Transformation matrix for mean to true equator and equinox
        l: Delaunay element (rad)
        l1: Delaunay element (rad)
        f: Delaunay element (rad)
        d: Delaunay element (rad)
        omega: Delaunay element (rad)
        lonmer: Longitude of Mercury (rad)
        lonven: Longitude of Venus (rad)
        lonear: Longitude of Earth (rad)
        lonmar: Longitude of Mars (rad)
        lonjup: Longitude of Jupiter (rad)
        lonsat: Longitude of Saturn (rad)
        lonurn: Longitude of Uranus (rad)
        lonnep: Longitude of Neptune (rad)
        precrate: Precession rate (rad/century)
    """
    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    deg2rad = np.pi / 180.0

    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2

    # Obtain data for calculations from the 2000B theory
    opt = '02'  # a-all, r-reduced, e-1980 theory
    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = fundarg(ttt, opt)

    # Obtain data coefficients
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()

    pnsum = 0.0
    ensum = 0.0
    for i in range(77, 0, -1):
        tempval = apni[i - 1, 0] * l + apni[i - 1, 1] * l1 + apni[i - 1, 2] * f + apni[i - 1, 3] * d + apni[i - 1, 4] * omega
        pnsum += (apn[i - 1, 0] + apn[i - 1, 1] * ttt) * np.sin(tempval) + (apn[i - 1, 4] + apn[i - 1, 5] * ttt) * np.cos(tempval)
        ensum += (apn[i - 1, 2] + apn[i - 1, 3] * ttt) * np.cos(tempval) + (apn[i - 1, 6] + apn[i - 1, 7] * ttt) * np.sin(tempval)

    # Form the planetary arguments
    pplnsum = -0.000135 * convrt  # " to rad
    eplnsum = 0.000388 * convrt

    # Add planetary and luni-solar components
    deltapsi = pnsum + pplnsum
    deltaeps = ensum + eplnsum

    prec, psia, wa, ea, xa = precess(ttt, '06')

    oblo = 84381.406 * convrt  # " to rad

    # Find nutation matrix
    a1 = rot1mat(ea + deltaeps)
    a2 = rot3mat(deltapsi)
    a3 = rot1mat(-ea)
    a4 = rot3mat(-xa)
    a5 = rot1mat(wa)
    a6 = rot3mat(psia)
    a7 = rot1mat(-oblo)
    a8 = rot1mat(-0.0068192 * convrt)
    a9 = rot2mat(0.0417750 * np.sin(oblo) * convrt)
    a10 = rot3mat(0.0146 * convrt)

    pnb = a10 @ a9 @ a8 @ a7 @ a6 @ a5 @ a4 @ a3 @ a2 @ a1
    nut = a3 @ a2 @ a1

    return deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate
