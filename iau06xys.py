import numpy as np
from fundarg import fundarg
from iau06in import iau06in

def iau06xys(ttt, ddx, ddy):
    """
    Calculates the transformation matrix that accounts for the
    effects of precession-nutation in the IAU2006 theory.

    Author: David Vallado, 719-573-2600, 16 Jul 2004

    Revisions:
    - Vallado: Consolidate with IAU 2000, 14 Feb 2005

    Inputs:
        ttt: Julian centuries of TT
        ddx: EOP correction for x (rad)
        ddy: EOP correction for y (rad)

    Outputs:
        x: Coordinate of CIP (rad)
        y: Coordinate of CIP (rad)
        s: Coordinate (rad)
        nut: Transformation matrix for TIRS-GCRF

    Locals:
        axs0: Real coefficients for x (rad)
        a0xi: Integer coefficients for x
        ays0: Real coefficients for y (rad)
        a0yi: Integer coefficients for y
        ass0: Real coefficients for s (rad)
        a0si: Integer coefficients for s
        apn: Real coefficients for nutation (rad)
        apni: Integer coefficients for nutation
        appl: Real coefficients for planetary nutation (rad)
        appli: Integer coefficients for planetary nutation
        ttt2, ttt3: Powers of ttt
        l, l1, f, d, omega: Delaunay elements (rad)
        deltaeps: Change in obliquity (rad)
        many others

    Coupling:
        iau00in: Initialize the arrays
        fundarg: Find the fundamental arguments

    References:
        Vallado 2004, 212-214
    """

    convrt = np.pi / (180.0 * 3600.0)
    deg2rad = np.pi / 180.0

    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2

    # Call iau06in function to initialize coefficients
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, _, _ = iau06in()

    opt = '06'  # 02 - 2000a, 96 - 1996 theory, 80-1980 theory

    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = fundarg( ttt, opt )
    x = -0.016617 + 2004.191898 * ttt - 0.4297829 * ttt2 - 0.19861834 * ttt3 - 0.000007578 * ttt4 + 0.0000059285 * ttt5
    y = -0.006951 - 0.025896 * ttt - 22.4072747 * ttt2 + 0.00190059 * ttt3 + 0.001112526 * ttt4 + 0.0000001358 * ttt5
    s = 0.000094 + 0.00380865 * ttt - 0.00012268 * ttt2 - 0.07257411 * ttt3 + 0.00002798 * ttt4 + 0.00001562 * ttt5

    # Initialize sums
    xsum0, xsum1, xsum2, xsum3, xsum4 = 0.0, 0.0, 0.0, 0.0, 0.0
    ysum0, ysum1, ysum2, ysum3, ysum4 = 0.0, 0.0, 0.0, 0.0, 0.0
    ssum0, ssum1, ssum2, ssum3, ssum4 = 0.0, 0.0, 0.0, 0.0, 0.0

    # Calculate xsum0
    for i in range(1306, 1, -1):
        tempval = (
            a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f + a0xi[i, 3] * d +
            a0xi[i, 4] * omega + a0xi[i, 5] * lonmer + a0xi[i, 6] * lonven +
            a0xi[i, 7] * lonear + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup +
            a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn + a0xi[i, 12] * lonnep +
            a0xi[i, 13] * precrate
        )
        xsum0 += axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    # Calculate xsum1
    for j in range(253, 1, -1):
        i = 1306 + j
        tempval = (
            a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f + a0xi[i, 3] * d +
            a0xi[i, 4] * omega + a0xi[i, 5] * lonmer + a0xi[i, 6] * lonven +
            a0xi[i, 7] * lonear + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup +
            a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn + a0xi[i, 12] * lonnep +
            a0xi[i, 13] * precrate
        )
        xsum1 += axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    # Calculate xsum2
    for j in range(36, 1, -1):
        i = 1306 + 253 + j
        tempval = (
            a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f + a0xi[i, 3] * d +
            a0xi[i, 4] * omega + a0xi[i, 5] * lonmer + a0xi[i, 6] * lonven +
            a0xi[i, 7] * lonear + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup +
            a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn + a0xi[i, 12] * lonnep +
            a0xi[i, 13] * precrate
        )
        xsum2 += axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    # Calculate xsum3
    for j in range(4, 1, -1):
        i = 1306 + 253 + 36 + j
        tempval = (
            a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f + a0xi[i, 3] * d +
            a0xi[i, 4] * omega + a0xi[i, 5] * lonmer + a0xi[i, 6] * lonven +
            a0xi[i, 7] * lonear + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup +
            a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn + a0xi[i, 12] * lonnep +
            a0xi[i, 13] * precrate
        )
        xsum3 += axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    # Calculate xsum4
    for j in range(1, 1, -1):
        i = 1306 + 253 + 36 + 4 + j
        tempval = (
            a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f + a0xi[i, 3] * d +
            a0xi[i, 4] * omega + a0xi[i, 5] * lonmer + a0xi[i, 6] * lonven +
            a0xi[i, 7] * lonear + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup +
            a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn + a0xi[i, 12] * lonnep +
            a0xi[i, 13] * precrate
        )
        xsum4 += axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    # Calculate ysum0
    for i in range(962, 1, -1):
        tempval = (
            a0yi[i, 0] * l + a0yi[i, 1] * l1 + a0yi[i, 2] * f + a0yi[i, 3] * d +
            a0yi[i, 4] * omega + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven +
            a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar + a0yi[i, 9] * lonjup +
            a0yi[i, 10] * lonsat + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep +
            a0yi[i, 13] * precrate
        )
        ysum0 += ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    # Calculate ysum1
    for j in range(277, 1, -1):
        i = 962 + j
        tempval = (
            a0yi[i, 0] * l + a0yi[i, 1] * l1 + a0yi[i, 2] * f + a0yi[i, 3] * d +
            a0yi[i, 4] * omega + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven +
            a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar + a0yi[i, 9] * lonjup +
            a0yi[i, 10] * lonsat + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep +
            a0yi[i, 13] * precrate
        )
        ysum1 += ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    # Calculate ysum2
    for j in range(30, 1, -1):
        i = 962 + 277 + j
        tempval = (
            a0yi[i, 0] * l + a0yi[i, 1] * l1 + a0yi[i, 2] * f + a0yi[i, 3] * d +
            a0yi[i, 4] * omega + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven +
            a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar + a0yi[i, 9] * lonjup +
            a0yi[i, 10] * lonsat + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep +
            a0yi[i, 13] * precrate
        )
        ysum2 += ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    # Calculate ysum3
    for j in range(5, 1, -1):
        i = 962 + 277 + 30 + j
        tempval = (
            a0yi[i, 0] * l + a0yi[i, 1] * l1 + a0yi[i, 2] * f + a0yi[i, 3] * d +
            a0yi[i, 4] * omega + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven +
            a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar + a0yi[i, 9] * lonjup +
            a0yi[i, 10] * lonsat + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep +
            a0yi[i, 13] * precrate
        )
        ysum3 += ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    # Calculate ysum4
    for j in range(1, 1, -1):
        i = 962 + 277 + 30 + 5 + j
        tempval = (
            a0yi[i, 0] * l + a0yi[i, 1] * l1 + a0yi[i, 2] * f + a0yi[i, 3] * d +
            a0yi[i, 4] * omega + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven +
            a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar + a0yi[i, 9] * lonjup +
            a0yi[i, 10] * lonsat + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep +
            a0yi[i, 13] * precrate
        )
        ysum4 += ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    # Calculate x and y
    x = x * convrt + xsum0 + xsum1 * ttt + xsum2 * ttt2 + xsum3 * ttt3 + xsum4 * ttt4
    y = y * convrt + ysum0 + ysum1 * ttt + ysum2 * ttt2 + ysum3 * ttt3 + ysum4 * ttt4

    for i in range(33, 1, -1):
        tempval = (
            a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2] * f + a0si[i, 3] * d +
            a0si[i, 4] * omega + a0si[i, 5] * lonmer + a0si[i, 6] * lonven +
            a0si[i, 7] * lonear + a0si[i, 8] * lonmar + a0si[i, 9] * lonjup +
            a0si[i, 10] * lonsat + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep +
            a0si[i, 13] * precrate
        )
        ssum0 += ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    # Calculate ssum1
    for j in range(3, 1, -1):
        i = 33 + j
        tempval = (
            a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2] * f + a0si[i, 3] * d +
            a0si[i, 4] * omega + a0si[i, 5] * lonmer + a0si[i, 6] * lonven +
            a0si[i, 7] * lonear + a0si[i, 8] * lonmar + a0si[i, 9] * lonjup +
            a0si[i, 10] * lonsat + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep +
            a0si[i, 13] * precrate
        )
        ssum1 += ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    # Calculate ssum2
    for j in range(25, 1, -1):
        i = 33 + 3 + j
        tempval = (
            a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2] * f + a0si[i, 3] * d +
            a0si[i, 4] * omega + a0si[i, 5] * lonmer + a0si[i, 6] * lonven +
            a0si[i, 7] * lonear + a0si[i, 8] * lonmar + a0si[i, 9] * lonjup +
            a0si[i, 10] * lonsat + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep +
            a0si[i, 13] * precrate
        )
        ssum2 += ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    # Calculate ssum3
    for j in range(4, 1, -1):
        i = 33 + 3 + 25 + j
        tempval = (
            a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2] * f + a0si[i, 3] * d +
            a0si[i, 4] * omega + a0si[i, 5] * lonmer + a0si[i, 6] * lonven +
            a0si[i, 7] * lonear + a0si[i, 8] * lonmar + a0si[i, 9] * lonjup +
            a0si[i, 10] * lonsat + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep +
            a0si[i, 13] * precrate
        )
        ssum3 += ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    # Calculate ssum4
    for j in range(1, 1, -1):
        i = 33 + 3 + 25 + 4 + j
        tempval = (
            a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2] * f + a0si[i, 3] * d +
            a0si[i, 4] * omega + a0si[i, 5] * lonmer + a0si[i, 6] * lonven +
            a0si[i, 7] * lonear + a0si[i, 8] * lonmar + a0si[i, 9] * lonjup +
            a0si[i, 10] * lonsat + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep +
            a0si[i, 13] * precrate
        )
        ssum4 += ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    # Calculate x, y, and s - all in radians
    x = x * convrt + xsum0 + xsum1 * ttt + xsum2 * ttt2 + xsum3 * ttt3 + xsum4 * ttt4
    y = y * convrt + ysum0 + ysum1 * ttt + ysum2 * ttt2 + ysum3 * ttt3 + ysum4 * ttt4
    s = -x * y * 0.5 + s * convrt + ssum0 + ssum1 * ttt + ssum2 * ttt2 + ssum3 * ttt3 + ssum4 * ttt4

    # Apply corrections
    x += ddx
    y += ddy

    # Now find a
    a = 0.5 + 0.125 * (x * x + y * y) # units take on whatever x and y are

    # Find nutation matrix
    nut1 = np.zeros((3, 3))
    nut1[0, 0] = 1.0 - a * x * x
    nut1[0, 1] = -a * x * y
    nut1[0, 2] = x
    nut1[1, 0] = -a * x * y
    nut1[1, 1] = 1.0 - a * y * y
    nut1[1, 2] = y
    nut1[2, 0] = -x
    nut1[2, 1] = -y
    nut1[2, 2] = 1.0 - a * (x * x + y * y)

    nut2 = np.eye(3)
    nut2[0, 0] = np.cos(s)
    nut2[1, 1] = np.cos(s)
    nut2[0, 1] = np.sin(s)
    nut2[1, 0] = -np.sin(s)

    nut = nut1.dot(nut2)

    return x, y, s, nut
