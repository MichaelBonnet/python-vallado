import math

def newtonnu(ecc, nu):
    """
    This function solves Kepler's equation when the true anomaly is known.
    The mean and eccentric, parabolic, or hyperbolic anomaly is also found.
    The parabolic limit at 168 is arbitrary. The hyperbolic anomaly is also
    limited. The hyperbolic sine is used because it's not double valued.

    Author: David Vallado, 719-573-2600, 27 May 2002

    Revisions:
    - Vallado: Fix small, 24 Sep 2002

    Inputs:
        - ecc: Eccentricity (0.0 to)
        - nu: True anomaly (-2pi to 2pi rad)

    Outputs:
        - e0: Eccentric anomaly (0.0 to 2pi rad, 153.02 deg)
        - m: Mean anomaly (0.0 to 2pi rad, 151.7425 deg)

    Locals:
        - e1: Eccentric anomaly, next value (rad)
        - sine: Sine of e
        - cose: Cosine of e
        - ktr: Index

    Coupling:
        - arcsinh: Arc hyperbolic sine
        - sinh: Hyperbolic sine

    References:
        - Vallado, 2007, 85, Algorithm 5
    """
    # Constants
    small = 0.00000001
    
    # Initialize variables
    e0 = 999999.9
    m = 999999.9

    # Circular case
    if abs(ecc) < small:
        m = nu
        e0 = nu
    else:
        # Elliptical case
        if ecc < 1.0 - small:
            sine = (math.sqrt(1.0 - ecc * ecc) * math.sin(nu)) / (1.0 + ecc * math.cos(nu))
            cose = (ecc + math.cos(nu)) / (1.0 + ecc * math.cos(nu))
            e0 = math.atan2(sine, cose)
            m = e0 - ecc * math.sin(e0)
        else:
            # Hyperbolic case
            if ecc > 1.0 + small:
                if ecc > 1.0 and abs(nu) + 0.00001 < math.pi - math.acos(1.0 / ecc):
                    sine = (math.sqrt(ecc * ecc - 1.0) * math.sin(nu)) / (1.0 + ecc * math.cos(nu))
                    e0 = math.asinh(sine)
                    m = ecc * math.sinh(e0) - e0
            else:
                # Parabolic case
                if abs(nu) < 168.0 * math.pi / 180.0:
                    e0 = math.tan(nu * 0.5)
                    m = e0 + (e0 * e0 * e0) / 3.0

    # Adjustments for eccentricity < 1
    if ecc < 1.0:
        m %= 2.0 * math.pi
        if m < 0.0:
            m += 2.0 * math.pi
        e0 %= 4.0 * math.pi

    return e0, m
