import numpy as np
from constastro import re, velkmps

def biellip(rinit, rb, rfinal, einit, efinal, nuinit, nufinal):
    """
    This function calculates the delta v's for a bi-elliptic transfer for either
    circle to circle, or ellipse to ellipse.

    Inputs:
        rinit : float
            Initial position magnitude (er)
        rb : float
            Interim orbit magnitude (er)
        rfinal : float
            Final position magnitude (er)
        einit : float
            Eccentricity of first orbit
        efinal : float
            Eccentricity of final orbit
        nuinit : float
            True anomaly of first orbit (0 or pi rad)
        nufinal : float
            True anomaly of final orbit (0 or pi rad, opposite of nuinit)

    Outputs:
        deltava : float
            Change in velocity at point a (er / tu)
        deltavb : float
            Change in velocity at point b (er / tu)
        deltavc : float
            Change in velocity at point c (er / tu)
        dttu : float
            Time of flight for the transfer (tu)

    """

    # --------------------  initialize values   ------------------- }
    mu = 1.0  # canonical units

    ainit = (rinit * (1.0 + einit * np.cos(nuinit))) / (1.0 - einit * einit)
    atran1 = (rinit + rb) * 0.5
    atran2 = (rb + rfinal) * 0.5
    afinal = (rfinal * (1.0 + efinal * np.cos(nufinal))) / (1.0 - efinal * efinal)

    deltava = 0.0
    deltavb = 0.0
    deltavc = 0.0
    dttu = 0.0

    if (einit < 1.0) and (efinal < 1.0):

        # -----------------  find delta v at point a  ----------------- }
        vinit = np.sqrt((2.0 * mu) / rinit - (mu / ainit))
        vtran1a = np.sqrt((2.0 * mu) / rinit - (mu / atran1))
        deltava = abs(vtran1a - vinit)

        # -----------------  find delta v at point b  ----------------- }
        vtran1b = np.sqrt((2.0 * mu) / rb - (mu / atran1))
        vtran2b = np.sqrt((2.0 * mu) / rb - (mu / atran2))
        deltavb = abs(vtran1b - vtran2b)

        # -----------------  find delta v at point c  ----------------- }
        vtran2c = np.sqrt((2.0 * mu) / rfinal - (mu / atran2))
        vfinal = np.sqrt((2.0 * mu) / rfinal - (mu / afinal))
        deltavc = abs(vfinal - vtran2c)

        # ----------------  find transfer time of flight  ------------- }
        dttu = np.pi * np.sqrt((atran1 * atran1 * atran1) / mu) + \
               np.pi * np.sqrt((atran2 * atran2 * atran2) / mu)

        # constastro;
        t1 = np.pi * np.sqrt((atran1 * atran1 * atran1) / 1) * 13.446852064
        t2 = np.pi * np.sqrt((atran2 * atran2 * atran2) / 1) * 13.446852064
        print(f' atran1   {atran1:11.7f}  {atran1 * re:11.7f} km')
        print(f' atran2   {atran2:11.7f}  {atran2 * re:11.7f} km')
        print(f' vinit   {vinit:11.7f}  {vinit * velkmps:11.7f} km/s')
        print(f' vtran1a  {vtran1a:11.7f}  {vtran1a * velkmps:11.7f} km/s')
        print(f' vtran1b  {vtran1b:11.7f}  {vtran1b * velkmps:11.7f} km/s')
        print(f' vtran2b  {vtran2b:11.7f}  {vtran2b * velkmps:11.7f} km/s')
        print(f' vtran2c  {vtran2c:11.7f}  {vtran2c * velkmps:11.7f} km/s')
        print(f' vfinal  {vfinal:11.7f}  {vfinal * velkmps:11.7f} km/s')
        print(f' t1  {t1:11.7f} t2 {t2:11.7f} min')

