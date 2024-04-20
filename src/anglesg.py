import numpy as np

def anglesg(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rs1, rs2, rs3):
    """
    This function solves the problem of orbit determination using three optical sightings.
    The solution function uses the Gaussian technique.

    Args:
        decl1, decl2, decl3: Declination for sightings 1, 2, and 3 (rad)
        rtasc1, rtasc2, rtasc3: Right ascension for sightings 1, 2, and 3 (rad)
        jd1, jdf1, jd2, jdf2, jd3, jdf3: Julian date of sightings
        rs1, rs2, rs3: ECI site position vectors (km)

    Returns:
        r2: IJK position vector at t2 (km)
        v2: IJK velocity vector at t2 (km/s)
    """
    # Constants
    re = 6378.137  # Radius of Earth (km)
    mu = 398600.4418  # Gravitational parameter of Earth (km^3/s^2)

    # Implementation
    ddpsi = 0.0  # delta corrections not needed for this level of precision
    ddeps = 0.0
    magr1in = 2.0 * re  # initial guesses
    magr2in = 2.01 * re
    direct = 'y'  # direction of motion short way

    # Set middle to 0, find decls to others
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0  # days to sec
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0

    # Find line of sight unit vectors
    los1 = np.array([np.cos(decl1) * np.cos(rtasc1),
                     np.cos(decl1) * np.sin(rtasc1),
                     np.sin(decl1)])

    los2 = np.array([np.cos(decl2) * np.cos(rtasc2),
                     np.cos(decl2) * np.sin(rtasc2),
                     np.sin(decl2)])

    los3 = np.array([np.cos(decl3) * np.cos(rtasc3),
                     np.cos(decl3) * np.sin(rtasc3),
                     np.sin(decl3)])

    # Find l matrix and determinant
    l1eci = los1
    l2eci = los2
    l3eci = los3

    lmat = np.vstack((l1eci, l2eci, l3eci)).T
    rsmat = np.vstack((rs1, rs2, rs3)).T

    li = np.linalg.inv(lmat)

    lcmat = np.zeros((3, 3))
    p1 = np.cross(los2, los3)
    p2 = np.cross(los1, los3)
    p3 = np.cross(los1, los2)

    dx = np.dot(los1, p1)
    lcmat[:, 0] = np.dot(rsmat, p1)
    lcmat[:, 1] = np.dot(rsmat, p2)
    lcmat[:, 2] = np.dot(rsmat, p3)

    tau31 = (jd3 - jd1) * 86400.0

    aa = 1 / dx * (-lcmat[0, 1] * tau32 / tau31 + lcmat[1, 1] + lcmat[2, 1] * tau12 / tau31)
    bb = 1 / (6.0 * dx) * (lcmat[0, 1] * (tau32 ** 2 - tau31 ** 2) * tau32 / tau31
                           + lcmat[2, 1] * (tau31 ** 2 - tau12 ** 2) * tau12 / tau31)

    d = np.linalg.det(lmat)
    lmati = np.linalg.inv(lmat)
    lir = np.dot(lmati, rsmat)

    # Form initial guess of r2
    dl1 = lir[1, 0] * aa - lir[1, 1] + lir[1, 2] * bb
    dl2 = lir[1, 0] * aa + lir[1, 2] * bb

    magrs2 = np.linalg.norm(rs2)
    l2dotrs = np.dot(los2, rs2)

    poly = np.zeros(9)
    poly[2] = -(dl1 ** 2 + 2.0 * dl1 * l2dotrs + magrs2 ** 2)
    poly[5] = -2.0 * mu * (l2dotrs * dl2 + dl1 * dl2)
    poly[8] = -mu ** 2 * dl2 ** 2

    rootarr = np.roots(poly)
    bigr2 = max([root.real for root in rootarr if root.imag == 0])

    # Solve matrix with u2 better known
    u = mu / bigr2 ** 3
    c1 = aa + aa * u
    c3 = bb + bb * u

    cmat = np.array([[-c1], [-1.0], [-c3]])
    rhomat = np.dot(lir, cmat)

    rhoold1 = rhomat[0, 0] / c1
    rhoold2 = rhomat[1, 0] / (-1.0)
    rhoold3 = rhomat[2, 0] / c3

    r1 = rhomat[0, 0] * l1eci / c1 + rs1
    r2 = rhomat[1, 0] * l2eci / (-1.0) + rs2
    r3 = rhomat[2, 0] * l3eci / c3 + rs3

    # Refine the answer
    rho2 = 999999.9
    ll = 0
    while abs(rhoold2 - rho2) > 1.0e-12 and ll <= 0:
        ll += 1
        rho2 = rhoold2

        r1 = rhomat[0, 0] * l1eci / c1 + rs1
        r2 = -rhomat[1, 0] * l2eci + rs2
        r3 = rhomat[2, 0] * l3eci / c3 + rs3

        magr1 = np.linalg.norm(r1)
        magr2 = np.linalg.norm(r2)
        magr3 = np.linalg.norm(r3)

        n1 = np.cross(r1, r2)
        n2 = np.cross(r2, r3)
        n3 = np.cross(r1, r3)

        mag_n1 = np.linalg.norm(n1)
        mag_n2 = np.linalg.norm(n2)
        mag_n3 = np.linalg.norm(n3)

        cosd12 = np.dot(r1, r2) / (magr1 * magr2)
        cosd23 = np.dot(r2, r3) / (magr2 * magr3)
        cosd31 = np.dot(r1, r3) / (magr1 * magr3)

        sind12 = mag_n1 / (magr1 * magr2)
        sind23 = mag_n2 / (magr2 * magr3)
        sind31 = mag_n3 / (magr1 * magr3)

        cos2d12 = np.dot(r1, r2) / (magr1 * magr2)
        cos2d23 = np.dot(r2, r3) / (magr2 * magr3)
        cos2d31 = np.dot(r1, r3) / (magr1 * magr3)

        s1sq = 1.0 - cosd12 ** 2
        s2sq = 1.0 - cosd23 ** 2
        s3sq = 1.0 - cosd31 ** 2

        u = np.zeros(3)
        u[0] = np.sqrt(mu) / magr1 ** 1.5
        u[1] = np.sqrt(mu) / magr2 ** 1.5
        u[2] = np.sqrt(mu) / magr3 ** 1.5

        f1 = -tau32 / tau31
        f2 = 1.0
        f3 = tau12 / tau31

        g1 = f2 * f3 / s1sq
        g2 = -f1 * f3 / s2sq
        g3 = f1 * f2 / s3sq

        rho1 = (s2sq * g1 * rhoold1 + s3sq * g2 * rhoold2 + s1sq * g3 * rhoold3) / \
               (s1sq * (g1 + g3) + s2sq * (g1 + g2) + s3sq * (g2 + g3))

        rho2 = (rhoold2 - (rho1 - rhoold1) * f1) / f2
        rho3 = (rhoold3 - (rho1 - rhoold1) * f3) / f3

        rhomat = np.array([[rho1], [-rho2], [rho3]])

    # Solve for velocity
    drdot = np.zeros(3)
    drdot[0] = (-rho1 * sind12) / (tau31 * (s1sq * s2sq) ** 0.5)
    drdot[1] = (rho1 * (cosd12 - cosd31 * cosd23)) / (tau31 * (s1sq * s2sq) ** 0.5)
    drdot[2] = (rho1 * sind23) / (tau31 * (s1sq * s2sq) ** 0.5)

    vv = np.linalg.inv(li) @ drdot

    # velocity components
    v2 = vv + rs2

    return r2, v2