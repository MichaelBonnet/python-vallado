import numpy as np
from mag import mag

def checkhitearth(altPad, r1, v1t, r2, v2t, nrev):
    """
    This function checks to see if the trajectory hits the earth during the transfer.

    Args:
    - altPad (float): pad for alt above surface (km)
    - r1 (list): initial position vector of int (km)
    - v1t (list): initial velocity vector of trns (km/s)
    - r2 (list): final position vector of int (km)
    - v2t (list): final velocity vector of trns (km/s)
    - nrev (int): number of revolutions (0, 1, 2, ...)

    Returns:
    - hitearth (str): 'y' if earth was impacted, 'n' otherwise
    - hitearthstr (str): "y - radii" if earth was impacted, "no" otherwise
    """
    
    # -------------------------- implementation -----------------
    show = 'n'
    mu = 3.986004418e5
    
    hitearth = 'n'
    hitearthstr = 'no'
    
    magr1 = mag(r1)
    magr2 = mag(r2)
    rpad = 6378.137 + altPad
    
    # check whether Lambert transfer trajectory hits the Earth
    # this check may not be needed depending on input data
    if (magr1 < rpad or magr2 < rpad):
        # hitting earth already at start or stop point
        hitearth = 'y'
        hitearthstr = hitearth + ' initradii'
        if show == 'y':
            print('hitearth?', hitearthstr)
    else:
        rdotv1 = np.dot(r1, v1t)
        rdotv2 = np.dot(r2, v2t)
        
        # Solve for a
        ainv = 2.0 / magr1 - mag(v1t)**2 / mu
        
        # Find ecos(E)
        ecosea1 = 1.0 - magr1 * ainv
        ecosea2 = 1.0 - magr2 * ainv
        
        # Determine radius of perigee
        # 4 distinct cases pass thru perigee
        # nrev > 0 you have to check
        if (nrev > 0):
            a = 1.0 / ainv
            # elliptical orbit
            if (a > 0.0):
                esinea1 = rdotv1 / np.sqrt(mu * a)
                ecc = np.sqrt(ecosea1 * ecosea1 + esinea1 * esinea1)
            else:
                # hyperbolic orbit
                esinea1 = rdotv1 / np.sqrt(mu * abs(-a))
                ecc = np.sqrt(ecosea1 * ecosea1 - esinea1 * esinea1)
            
            rp = a * (1.0 - ecc)
            if (rp < rpad):
                hitearth = 'y'
                hitearthstr = hitearth + 'Sub_Earth_nrev'
        
        # nrev = 0, 3 cases:
        # heading to perigee and ending after perigee
        # both headed away from perigee, but end is closer to perigee
        # both headed toward perigee, but start is closer to perigee
        else:
            if ((rdotv1 < 0.0 and rdotv2 > 0.0) or (rdotv1 > 0.0 and rdotv2 > 0.0 and ecosea1 < ecosea2) or \
                    (rdotv1 < 0.0 and rdotv2 < 0.0 and ecosea1 > ecosea2)):
                # parabola
                if (abs(ainv) <= 1.0e-10):
                    hbar = np.cross(r1, v1t)
                    magh = mag(hbar)  # find h magnitude
                    rp = magh * magh * 0.5 / mu
                    if (rp < rpad):
                        hitearth = 'y'
                        hitearthstr = hitearth + ' Sub_Earth_para'
                else:
                    # for both elliptical & hyperbolic
                    a = 1.0 / ainv
                    esinea1 = rdotv1 / np.sqrt(mu * abs(a))
                    if (ainv > 0.0):
                        ecc = np.sqrt(ecosea1 * ecosea1 + esinea1 * esinea1)
                    else:
                        ecc = np.sqrt(ecosea1 * ecosea1 - esinea1 * esinea1)
                    if (ecc < 1.0):
                        rp = a * (1.0 - ecc)
                        if (rp < rpad):
                            hitearth = 'y'
                            hitearthstr = hitearth + ' Sub_Earth_ell'
                    else:
                        # hyperbolic heading towards the earth
                        if (rdotv1 < 0.0 and rdotv2 > 0.0):
                            rp = a * (1.0 - ecc)
                            if (rp < rpad):
                                hitearth = 'y'
                                hitearthstr = hitearth + ' Sub_Earth_hyp'
        
        if show == 'y':
            print('hitearth?', hitearth, 'rp', rp * 6378.0, 'km')

    return hitearth, hitearthstr
