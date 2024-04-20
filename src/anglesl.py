import numpy as np
from mag import mag
from constastro import earthrot, mu

# Function to solve the problem of orbit determination using three optical sightings and the method of Laplace
def anglesl(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rs1, rs2, rs3):
    """
    This function solves the problem of orbit determination using three optical sightings and the method of Laplace.

    Inputs:
        decl1, decl2, decl3: Declination #1, Declination #2, Declination #3 (rad)
        rtasc1, rtasc2, rtasc3: Right ascension #1, Right ascension #2, Right ascension #3 (rad)
        jd1, jdf1, jd2, jdf2, jd3, jdf3: Julian date of 1st, 2nd, and 3rd sightings (days from 4713 BC)
        diffsites: Flag indicating if sightings are from different sites ('n' for no, 'y' for yes)
        rs1, rs2, rs3: ECI site position vectors (km)

    Outputs:
        r2: IJK position vector (km)
        v2: IJK velocity vector (km/s)
    """

    # Constants
    earthrate = np.array([0.0, 0.0, earthrot])  # Earth's rotation rate (km/s)
    small = 0.00000001  # Tolerance
    rad = 180.0 / np.pi  # Conversion factor from radians to degrees

    # Time differences
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0  # Days to seconds
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0

    # Line of sight vectors
    los1 = np.array([np.cos(decl1) * np.cos(rtasc1), np.cos(decl1) * np.sin(rtasc1), np.sin(decl1)])
    los2 = np.array([np.cos(decl2) * np.cos(rtasc2), np.cos(decl2) * np.sin(rtasc2), np.sin(decl2)])
    los3 = np.array([np.cos(decl3) * np.cos(rtasc3), np.cos(decl3) * np.sin(rtasc3), np.sin(decl3)])

    # Lagrange interpolation coefficients
    s1 = -tau32 / (tau12 * tau13)
    s2 = (tau12 + tau32) / (tau12 * tau32)
    s3 = -tau12 / (-tau13 * tau32)
    s4 = 2.0 / (tau12 * tau13)
    s5 = 2.0 / (tau12 * tau32)
    s6 = 2.0 / (-tau13 * tau32)

    # Derivatives of line of sight
    ldot = s1 * los1 + s2 * los2 + s3 * los3  # rad / s
    lddot = s4 * los1 + s5 * los2 + s6 * los3  # rad / s^2

    # Second derivative of rs
    if diffsites == 'n':
        rs2dot = np.cross(earthrate, rs2)
        rs2ddot = np.cross(earthrate, rs2dot)
    else:
        rs2dot = s1 * rs1 + s2 * rs2 + s3 * rs3
        rs2ddot = s4 * rs1 + s5 * rs2 + s6 * rs3

    # Position and velocity determinant matrices
    dmat = np.vstack((los2, ldot, lddot))
    dmat1 = np.vstack((los2, ldot, rs2ddot))
    dmat2 = np.vstack((los2, ldot, rs2))
    dmat3 = np.vstack((los2, rs2ddot, lddot))
    dmat4 = np.vstack((los2, rs2, lddot))

    # Determinants
    d = 2.0 * np.linalg.det(dmat)
    d1 = np.linalg.det(dmat1)
    d2 = np.linalg.det(dmat2)
    d3 = np.linalg.det(dmat3)
    d4 = np.linalg.det(dmat4)

    # Determine rho magnitude
    if np.abs(d) > 1.0 - 1e-14:
        l2dotrs = np.dot(los2, rs2)
        poly = np.zeros(10)
        poly[2] = (l2dotrs * 4.0 * d1 / d - 4.0 * d1 * d1 / (d * d) - mag(rs2) ** 2)
        poly[5] = mu * (l2dotrs * 4.0 * d2 / d - 8.0 * d1 * d2 / (d * d))
        poly[8] = -4.0 * mu * mu * d2 * d2 / (d * d)
        rootarr = np.roots(poly)

        # Find the real root
        bigr2 = 0.0
        for root in rootarr:
            if root.real > bigr2:
                bigr2 = root.real

        rho = -2.0 * d1 / d - 2.0 * mu * d2 / (bigr2 ** 3 * d)

        # Find middle position vector
        r2 = rho * los2 + rs2

        # Find rhodot magnitude
        magr2 = mag(r2)
        rhodot = -d3 / d - mu * d4 / (magr2 ** 3 * d)

        # Find middle velocity vector
        v2 = rhodot * los2 + rho * ldot + rs2dot
    else:
        print("Determinant value was zero:", d)

    return r2, v2
