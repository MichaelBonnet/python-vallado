import numpy as np
from fundarg import fundarg
from iau06in import iau06in
from precess import precess
from rotations import rot1mat, rot2mat, rot3mat

def iau06pna(ttt):
    """
    Calculates the transformation matrix that accounts for the effects of precession-nutation in the IAU2000A theory.

    Parameters:
    - ttt : float
        Julian centuries of TT.

    Returns:
    - deltapsi : float
        Change in longitude (rad).
    - pnb : numpy.ndarray
        Nutation transformation matrix for IRE-GCRF.
    - prec : numpy.ndarray
        Precession transformation matrix.
    - nut : numpy.ndarray
        Nutation transformation matrix.
    - l : float
        Delaunay element (rad).
    - l1 : float
        Delaunay element (rad).
    - f : float
        Delaunay element (rad).
    - d : float
        Delaunay element (rad).
    - omega : float
        Delaunay element (rad).
    - lonmer : float
        Planetary value (rad).
    - lonven : float
        Planetary value (rad).
    - lonear : float
        Planetary value (rad).
    - lonmar : float
        Planetary value (rad).
    - lonjup : float
        Planetary value (rad).
    - lonsat : float
        Planetary value (rad).
    - lonurn : float
        Planetary value (rad).
    - lonnep : float
        Planetary value (rad).
    - precrate : float
        Planetary value (rad).
    """

    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    deg2rad = np.pi / 180.0

    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2

    # Obtain data for calculations from the 2000a theory
    opt = '06'  # a-all, r-reduced, e-1980 theory
    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = fundarg(ttt, opt)

    # ---- obtain data coefficients
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()

    pnsum = 0.0
    ensum = 0.0
    for i in range(678, 0, -1):
        tempval = apni[i - 1, 0] * l + apni[i - 1, 1] * l1 + apni[i - 1, 2] * f + apni[i - 1, 3] * d + apni[i - 1, 4] * omega
        tempval = np.mod(tempval, 2 * np.pi)  # rad
        pnsum += (apn[i - 1, 0] + apn[i - 1, 1] * ttt) * np.sin(tempval) + (apn[i - 1, 4]) * np.cos(tempval)
        ensum += (apn[i - 1, 2] + apn[i - 1, 3] * ttt) * np.cos(tempval) + (apn[i - 1, 6]) * np.sin(tempval)

    pplnsum = 0.0
    eplnsum = 0.0
    for i in range(1, 688):
        tempval = appli[i - 1, 0] * l + appli[i - 1, 1] * l1 + appli[i - 1, 2] * f + appli[i - 1, 3] * d + appli[i - 1, 4] * omega + \
                  appli[i - 1, 5] * lonmer + appli[i - 1, 6] * lonven + appli[i - 1, 7] * lonear + appli[i - 1, 8] * lonmar + \
                  appli[i - 1, 9] * lonjup + appli[i - 1, 10] * lonsat + appli[i - 1, 11] * lonurn + appli[i - 1, 12] * lonnep + \
                  appli[i - 1, 13] * precrate
        pplnsum += appl[i - 1, 0] * np.sin(tempval) + appl[i - 1, 1] * np.cos(tempval)
        eplnsum += appl[i - 1, 2] * np.sin(tempval) + appl[i - 1, 3] * np.cos(tempval)

    # Add planetary and luni-solar components
    deltapsi = pnsum + pplnsum  # rad
    deltaeps = ensum + eplnsum

    # Corrections to the IAU2000A
    j2d = -2.7774e-6 * ttt * convrt  # rad
    deltapsi += deltapsi * (0.4697e-6 + j2d)  # rad
    deltaeps += deltaeps * j2d

    prec, psia, wa, ea, xa = precess(ttt, '06')

    oblo = 84381.406 * convrt  # " to rad

    # Mean to true
    a1 = rot1mat(ea + deltaeps)
    a2 = rot3mat(deltapsi)
    a3 = rot1mat(-ea)

    # J2000 to date (precession)
    a4 = rot3mat(-xa)
    a5 = rot1mat(wa)
    a6 = rot3mat(psia)
    a7 = rot1mat(-oblo)

    # ICRS to J2000
    a8 = rot1mat(-0.0068192 * convrt)
    a9 = rot2mat(0.0417750 * np.sin(oblo) * convrt)
    a10 = rot3mat(0.0146 * convrt)

    pnb = a10.dot(a9).dot(a8).dot(a7).dot(a6).dot(a5).dot(a4).dot(a3).dot(a2).dot(a1)

    prec = a7.dot(a6).dot(a5).dot(a4)

    nut = a3.dot(a2).dot(a1)

    frb = a10.dot(a9).dot(a8)

    return deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
           lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate
