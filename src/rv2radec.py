import numpy as np

def rv2radec(r, v):
    """
    This function converts the right ascension and declination values with
    position and velocity vectors of a satellite. It uses the velocity vector to
    find the solution of singular cases.

    Inputs:
        r: Position vector (km)
        v: Velocity vector (km/s)

    Outputs:
        rr: Radius of the satellite (km)
        rtasc: Right ascension (rad)
        decl: Declination (rad)
        drr: Radius of the satellite rate (km/s)
        drtasc: Right ascension rate (rad/s)
        ddecl: Declination rate (rad/s)
    """
    small = 0.00000001
    
    # ------------- calculate angles and rates ----------------
    rr = np.linalg.norm(r)
    temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
    if temp < small:
        rtasc = np.arctan2(v[1], v[0])
    else:
        rtasc = np.arctan2(r[1], r[0])
    if rtasc < 0.0:
        rtasc += 2.0*np.pi
    decl = np.arcsin(r[2] / rr)

    temp1 = -r[1]*r[1] - r[0]*r[0]
    drr = np.dot(r, v) / rr
    if abs(temp1) > small:
        drtasc = (v[0]*r[1] - v[1]*r[0]) / temp1
    else:
        drtasc = 0.0
    if abs(temp) > small:
        ddecl = (v[2] - drr*np.sin(decl)) / temp
    else:
        ddecl = 0.0
    
    return rr, rtasc, decl, drr, drtasc, ddecl