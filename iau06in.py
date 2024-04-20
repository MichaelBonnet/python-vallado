import math

def iau06in():
    """
    Initialize matrices for IAU 2006 reduction calculations.

    Author: David Vallado, 719-573-2600, 16 Jul 2004
    Revisions:
    - Dav 14 Apr 11: Update for IAU2006 conventions

    Outputs:
        axs0 - Real coefficients for x (rad)
        a0xi - Integer coefficients for x
        ays0 - Real coefficients for y (rad)
        a0yi - Integer coefficients for y
        ass0 - Real coefficients for s (rad)
        a0si - Integer coefficients for s
        apn - Real coefficients for nutation (rad)
        apni - Integer coefficients for nutation
        appl - Real coefficients for obliquity (rad)
        appli - Integer coefficients for obliquity
        agst - Real coefficients for GST (rad)
        agsti - Integer coefficients for GST
    """

    # Conversion factors to radians
    convrtu = (0.000001 * math.pi) / (180.0 * 3600.0)  # If micro arcsecond
    convrtm = (0.001 * math.pi) / (180.0 * 3600.0)      # If milli arcsecond

    # XYS values
    filein = load('iau06xtab5.2.a.dat')
    axs0 = filein[:, 2:3]  # Reals
    a0xi = filein[:, 4:17]  # Integers
    for i in range(axs0.size):
        axs0[i, 0] *= convrtu  # rad
        axs0[i, 1] *= convrtu  # rad

    filein = load('iau06ytab5.2.b.dat')
    ays0 = filein[:, 2:3]
    a0yi = filein[:, 4:17]
    for i in range(ays0.size):
        ays0[i, 0] *= convrtu
        ays0[i, 1] *= convrtu

    filein = load('iau06stab5.2.d.dat')
    ass0 = filein[:, 2:3]
    a0si = filein[:, 4:17]
    for i in range(ass0.size):
        ass0[i, 0] *= convrtu
        ass0[i, 1] *= convrtu

    # Nutation values old approach IAU2003
    filein = load('iau03n.dat')
    apni = filein[:, 0:5]
    apn = filein[:, 6:13]
    for i in range(apn.size):
        apn[i, 0] *= convrtm
        apn[i, 1] *= convrtm
        apn[i, 2] *= convrtm
        apn[i, 3] *= convrtm
        apn[i, 4] *= convrtm
        apn[i, 5] *= convrtm
        apn[i, 6] *= convrtm
        apn[i, 7] *= convrtm

    # Planetary nutation values
    filein = load('iau03pl.dat')
    appli = filein[:, 1:14]
    appl = filein[:, 16:20]  # 21 is extra
    for i in range(appl.size):
        appl[i, 0] *= convrtm
        appl[i, 1] *= convrtm
        appl[i, 2] *= convrtm
        appl[i, 3] *= convrtm

    # GMST values
    filein = load('iau06gsttab5.2.e.dat')
    agst = filein[:, 1:3]
    agsti = filein[:, 3:16]
    for i in range(agst.size):
        agst[i, 0] *= convrtu
        agst[i, 1] *= convrtu

    return axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti
