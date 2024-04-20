import numpy as np
from mag import mag
from newtonnu import newtonnu

def rv2coe(r, v, mu):
    """
    Finds the classical orbital elements given the geocentric
    equatorial position and velocity vectors.

    Args:
        r (numpy.ndarray): Position vector in km.
        v (numpy.ndarray): Velocity vector in km/s.
        mu (float): Gravitational parameter in km^3/s^2.

    Returns:
        tuple: A tuple containing the following classical orbital elements:
            - p (float): Semilatus rectum in km.
            - a (float): Semimajor axis in km.
            - ecc (float): Eccentricity.
            - incl (float): Inclination in radians (0.0 to pi).
            - raan (float): Longitude of ascending node in radians (0.0 to 2pi).
            - argp (float): Argument of perigee in radians (0.0 to 2pi).
            - nu (float): True anomaly in radians (0.0 to 2pi).
            - m (float): Mean anomaly in radians (0.0 to 2pi).
            - arglat (float): Argument of latitude in radians (0.0 to 2pi).
            - truelon (float): True longitude in radians (0.0 to 2pi).
            - lonper (float): Longitude of periapsis in radians (0.0 to 2pi).
    """
    constmath = None  # Define your constants here
    constastro = None  # Define your constants here
    muin = mu  # this is the km version
    small = 1e-12
    
    # -------------------------  implementation   -----------------
    magr = mag(r)
    magv = mag(v)
    
    # ------------------  find h n and e vectors   ----------------
    hbar = np.cross(r, v)
    magh = mag(hbar)
    
    if magh >= 0.0:
        nbar = np.array([-hbar[1], hbar[0], 0.0])
        magn = mag(nbar)
        c1 = magv * magv - muin / magr
        rdotv = np.dot(r, v)
        
        ebar = (c1 * r - rdotv * v) / muin
        ecc = mag(ebar)
        
        # ------------  find a e and semi-latus rectum   ----------
        sme = magv * magv * 0.5 - muin / magr
        
        if abs(sme) > small:
            a = -muin / (2.0 * sme)
        else:
            a = np.inf
        
        p = magh * magh / muin
        
        # -----------------  find inclination   -------------------
        hk = hbar[2] / magh
        incl = np.arccos(hk)
        
        # --------  determine type of orbit for later use  --------
        typeorbit = 'ei'
        if ecc < small:
            if incl < small or abs(incl - np.pi) < small:
                typeorbit = 'ce'
            else:
                typeorbit = 'ci'
        else:
            if incl < small or abs(incl - np.pi) < small:
                typeorbit = 'ee'
        
        # ----------  find right ascension of ascending node ------------
        if magn > small:
            temp = nbar[0] / magn
            if abs(temp) > 1.0:
                temp = np.sign(temp)
            raan = np.arccos(temp)
            if nbar[1] < 0.0:
                raan = 2.0 * np.pi - raan
        else:
            raan = np.nan
        
        # ---------------- find argument of perigee ---------------
        if typeorbit == 'ei':
            argp = np.arccos(np.dot(nbar, ebar) / (mag(nbar) * mag(ebar)))
            if ebar[2] < 0.0:
                argp = 2.0 * np.pi - argp
        else:
            argp = np.nan
        
        # ------------  find true anomaly at epoch    -------------
        if typeorbit[0] == 'e':
            nu = np.arccos(np.dot(ebar, r) / (mag(ebar) * mag(r)))
            if rdotv < 0.0:
                nu = 2.0 * np.pi - nu
        else:
            nu = np.nan
        
        # ----  find argument of latitude - circular inclined -----
        # -- find in general cases too
        if typeorbit == 'ci' or typeorbit == 'ei':
            arglat = np.arccos(np.dot(nbar, r) / (mag(nbar) * mag(r)))
            if r[2] < 0.0:
                arglat = 2.0 * np.pi - arglat
            m = arglat
        else:
            arglat = np.nan
        
        # -- find longitude of perigee - elliptical equatorial ----
        if ecc > small and typeorbit == 'ee':
            temp = ebar[0] / ecc
            if abs(temp) > 1.0:
                temp = np.sign(temp)
            lonper = np.arccos(temp)
            if ebar[1] < 0.0:
                lonper = 2.0 * np.pi - lonper
            if incl > np.pi / 2.0:
                lonper = 2.0 * np.pi - lonper
        else:
            lonper = np.nan
        
        # -------- find true longitude - circular equatorial ------
        if magr > small and typeorbit == 'ce':
            temp = r[0] / magr
            if abs(temp) > 1.0:
                temp = np.sign(temp)
            truelon = np.arccos(temp)
            if r[1] < 0.0:
                truelon = 2.0 * np.pi - truelon
            if incl > np.pi / 2.0:
                truelon = 2.0 * np.pi - truelon
            m = truelon
        else:
            truelon = np.nan
        
        # ------------ find mean anomaly for all orbits -----------
        # if typeorbit[0] == 'e':
        e, m = newtonnu(ecc, nu)
        # end

    else:
        p = np.nan
        a = np.nan
        ecc = np.nan
        incl = np.nan
        raan = np.nan
        argp = np.nan
        nu = np.nan
        m = np.nan
        arglat = np.nan
        truelon = np.nan
        lonper = np.nan

    return p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper
