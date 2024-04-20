import numpy as np

def adbar2rv(rmag, vmag, rtasc, decl, fpav, az):
    """
    Function: adbar2rv

    This function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
    position and velocity magnitude) into ECI position and velocity vectors.

    Author: David Vallado (719-573-2600), 9 Jun 2002

    Inputs:
        rmag: ECI position vector magnitude (km)
        vmag: ECI velocity vector magnitude (km/sec)
        rtasc: Right ascension of satellite (rad)
        decl: Declination of satellite (rad)
        fpav: Satellite flight path angle from vertical (rad)
        az: Satellite flight path azimuth (rad)

    Outputs:
        r: ECI position vector (km)
        v: ECI velocity vector (km/s)

    Locals:
        None

    Coupling:
        None

    References:
        Vallado, 2001, xx
        Chobotov, 70

    Usage:
        r, v = adbar2rv(rmag, vmag, rtasc, decl, fpav, az)
    """


    # Form position vector
    r = np.zeros(3)
    r[0] = rmag * np.cos(decl) * np.cos(rtasc)
    r[1] = rmag * np.cos(decl) * np.sin(rtasc)
    r[2] = rmag * np.sin(decl)

    # Form velocity vector
    v = np.zeros(3)
    v[0] = vmag * (np.cos(rtasc) * (-np.cos(az) * np.sin(fpav) * np.sin(decl) + np.cos(fpav) * np.cos(decl)) - np.sin(az) * np.sin(fpav) * np.sin(rtasc))
    v[1] = vmag * (np.sin(rtasc) * (-np.cos(az) * np.sin(fpav) * np.sin(decl) + np.cos(fpav) * np.cos(decl)) + np.sin(az) * np.sin(fpav) * np.cos(rtasc))
    v[2] = vmag * (np.cos(az) * np.cos(decl) * np.sin(fpav) + np.cos(fpav) * np.sin(decl))

    return r, v
