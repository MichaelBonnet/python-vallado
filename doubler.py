import numpy as np

def doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
            los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct, re, mu):
    """
    This routine accomplishes the iteration work for the double-r angles only routine.

    Args:
    - cc1, cc2: constants
    - magrsite1, magrsite2: magnitudes
    - magr1in, magr2in: magnitudes
    - los1, los2, los3: vectors
    - rsite1, rsite2, rsite3: position vectors
    - t1, t3: time values
    - direct: 'y' or 'n'
    - re: Earth's radius
    - mu: gravitational parameter

    Returns:
    - r2, r3: position vectors
    - f1, f2: time values
    - q1: value
    - magr1, magr2, a, deltae32: magnitudes and angle
    """

    rho1 = (-cc1 + np.sqrt(cc1**2 - 4.0 * (magrsite1**2 - magr1in**2))) / 2.0
    rho2 = (-cc2 + np.sqrt(cc2**2 - 4.0 * (magrsite2**2 - magr2in**2))) / 2.0
    r1 = rho1 * los1 + rsite1
    r2 = rho2 * los2 + rsite2

    if direct == 'y':
        w = np.cross(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))
    else:
        w = -np.cross(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))

    rho3 = -np.dot(rsite3, w) / np.dot(los3, w)
    r3 = rho3 * los3 + rsite3
    magr1 = np.linalg.norm(r1)
    magr2 = np.linalg.norm(r2)
    magr3 = np.linalg.norm(r3)

    cosdv21 = np.dot(r2, r1) / (magr2 * magr1)
    sindv21 = np.linalg.norm(np.cross(r2, r1)) / (magr2 * magr1)
    dv21 = np.arctan2(sindv21, cosdv21)

    cosdv31 = np.dot(r3, r1) / (magr3 * magr1)
    sindv31 = np.sqrt(1.0 - cosdv31**2)
    dv31 = np.arctan2(sindv31, cosdv31)

    cosdv32 = np.dot(r3, r2) / (magr3 * magr2)
    sindv32 = np.linalg.norm(np.cross(r3, r2)) / (magr3 * magr2)
    dv32 = np.arctan2(sindv32, cosdv32)

    if dv31 > np.pi:
        c1 = (magr2 * sindv32) / (magr1 * sindv31)
        c3 = (magr2 * sindv21) / (magr3 * sindv31)
        p = (c1 * magr1 + c3 * magr3 - magr2) / (c1 + c3 - 1)
    else:
        c1 = (magr1 * sindv31) / (magr2 * sindv32)
        c3 = (magr1 * sindv21) / (magr3 * sindv32)
        p = (c3 * magr3 - c1 * magr2 + magr1) / (-c1 + c3 + 1)

    ecosv1 = p / magr1 - 1
    ecosv2 = p / magr2 - 1
    ecosv3 = p / magr3 - 1

    if dv21 != np.pi:
        esinv2 = (-cosdv21 * ecosv2 + ecosv1) / sindv21
    else:
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv32

    e = np.sqrt(ecosv2**2 + esinv2**2)
    a = p / (1 - e**2)

    if e**2 < 0.99:
        n = np.sqrt(mu / a**3)

        s = magr2 / p * np.sqrt(1 - e**2) * esinv2
        c = magr2 / p * (e**2 + ecosv2)

        sinde32 = magr3 / np.sqrt(a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s
        cosde32 = 1 - magr2 * magr3 / (a * p) * (1 - cosdv32)
        deltae32 = np.arctan2(sinde32, cosde32)

        sinde21 = magr1 / np.sqrt(a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s
        cosde21 = 1 - magr2 * magr1 / (a * p) * (1 - cosdv21)
        deltae21 = np.arctan2(sinde21, cosde21)

        deltam32 = deltae32 + 2 * s * (np.sin(deltae32 / 2))**2 - c * np.sin(deltae32)
        deltam12 = -deltae21 + 2 * s * (np.sin(deltae21 / 2))**2 + c * np.sin(deltae21)
    else:
        print('Hyperbolic, e1 is greater than 0.99:', e)
        n = np.sqrt(mu / (-a**3))

        s = magr2 / p * np.sqrt(e**2 - 1) * esinv2
        c = magr2 / p * (e**2 + ecosv2)

        sindh32 = magr3 / np.sqrt(-a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s
        deltah32 = np.log(sindh32 + np.sqrt(sindh32**2 + 1))

        sindh21 = magr1 / np.sqrt(-a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s
        deltah21 = np.log(sindh21 + np.sqrt(sindh21**2 + 1))

        deltam32 = -deltah32 + 2 * s * (np.sinh(deltah32 / 2))**2 + c * np.sinh(deltah32)
        deltam12 = deltah21 + 2 * s * (np.sinh(deltah21 / 2))**2 - c * np.sinh(deltah21)
        deltae32 = deltah32  # fix

    f1 = t1 - deltam12 / n
    f2 = t3 - deltam32 / n
    q1 = np.sqrt(f1**2 + f2**2)

    return r2, r3, f1, f2, q1, magr1, magr2, a, deltae32
