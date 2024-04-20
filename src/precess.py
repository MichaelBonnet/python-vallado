import numpy as np

def precess(ttt, opt):
    """
    This function calculates the transformation matrix that accounts for the effects
    of precession. Both the 1980 and 2006 theories are handled. Note that the
    required parameters differ a little.

    Author: David Vallado, 719-573-2600, 25 Jun 2002

    Revisions:
    - Vallado: Consolidate with iau 2000, 14 Feb 2005

    Inputs:
    - ttt: Julian centuries of TT
    - opt: Method option: '01', '02', '96', '80'

    Outputs:
    - prec: Transformation matrix for mod - j2000 (80 only)
    - psia: Canonical precession angle in radians (00 only)
    - wa: Canonical precession angle in radians (00 only)
    - ea: Canonical precession angle in radians (00 only)
    - xa: Canonical precession angle in radians (00 only)
    """
    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt

    prec = np.eye(3)

    # ------------------- fk4 b1950 precession angles --------------------
    if opt == '50':
        # Newcomb Exp Supp 61 approach
        t1 = 0.0
        t2 = (ttt - 2433282.42345905) / 36525  # uses B1950
        psia = 50.3708 + 0.0050 * ttt
        wa = 0.0  # not sure which one is which...
        ea = 84428.26 - 46.845 * ttt - 0.00059 * ttt2 + 0.00181 * ttt3
        xa = 0.1247 - 0.0188 * ttt
        # GTDS pg 3-17 using days from 1950 - avoids long precession
        zeta = 2304.9969 * ttt + 0.302 * ttt2 + 0.01808 * ttt3
        theta = 2004.2980 * ttt - 0.425936 * ttt2 - 0.0416 * ttt3
        z = 2304.9969 * ttt + 1.092999 * ttt2 + 0.0192 * ttt3
        # tp-008 36-45
        # ttt is tropical centuries from 1950 36524.22 days
        prec[0, 0] = 1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3
        prec[0, 1] = 2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3
        prec[0, 2] = 9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3
        prec[1, 0] = -prec[0, 1]
        prec[1, 1] = 1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3
        prec[1, 2] = -1.0858e-4 * ttt2
        prec[2, 0] = -prec[0, 2]
        prec[2, 1] = prec[1, 2]
        prec[2, 2] = 1.0 - 4.721e-5 * ttt2
        # pass these back out for testing
        psia = zeta
        wa = theta
        ea = z
    else:
        if opt == '80':
            psia = 5038.7784 * ttt - 1.07259 * ttt2 - 0.001147 * ttt3
            wa = 84381.448 + 0.05127 * ttt2 - 0.007726 * ttt3
            ea = 84381.448 - 46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3
            xa = 10.5526 * ttt - 2.38064 * ttt2 - 0.001125 * ttt3
            zeta = 2306.2181 * ttt + 0.30188 * ttt2 + 0.017998 * ttt3
            theta = 2004.3109 * ttt - 0.42665 * ttt2 - 0.041833 * ttt3
            z = 2306.2181 * ttt + 1.09468 * ttt2 + 0.018203 * ttt3
        else:
            oblo = 84381.406
            psia = (((( -0.0000000951 * ttt + 0.000132851 ) * ttt - 0.00114045 ) * ttt - 1.0790069 ) * ttt + 5038.481507 ) * ttt
            wa = ((((  0.0000003337 * ttt - 0.000000467 ) * ttt - 0.00772503 ) * ttt + 0.0512623 ) * ttt - 0.025754 ) * ttt + oblo
            ea = (((( -0.0000000434 * ttt - 0.000000576 ) * ttt + 0.00200340 ) * ttt - 0.0001831 ) * ttt - 46.836769 ) * ttt + oblo
            xa = (((( -0.0000000560 * ttt + 0.000170663 ) * ttt - 0.00121197 ) * ttt - 2.3814292 ) * ttt + 10.556403 ) * ttt
            zeta = (((( -0.0000003173 * ttt - 0.000005971 ) * ttt + 0.01801828 ) * ttt + 0.2988499 ) * ttt + 2306.083227 ) * ttt + 2.650545
            theta = (((( -0.0000001274 * ttt - 0.000007089 ) * ttt - 0.04182264 ) * ttt - 0.4294934 ) * ttt + 2004.191903 ) * ttt
            z = ((((  0.0000002904 * ttt - 0.000028596 ) * ttt + 0.01826837 ) * ttt + 1.0927348 ) * ttt + 2306.077181 ) * ttt - 2.650545

    # convert units to rad
    psia = psia * convrt  # rad
    wa = wa * convrt
    ea = ea * convrt
    xa = xa * convrt

    zeta = zeta * convrt
    theta = theta * convrt
    z = z * convrt

    if opt == '80' or opt == '06':
        coszeta = np.cos(zeta)
        sinzeta = np.sin(zeta)
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        cosz = np.cos(z)
        sinz = np.sin(z)
        # ----------------- form matrix  mod to j2000 -----------------
        prec[0, 0] = coszeta * costheta * cosz - sinzeta * sinz
        prec[0, 1] = coszeta * costheta * sinz + sinzeta * cosz
        prec[0, 2] = coszeta * sintheta
        prec[1, 0] = -sinzeta * costheta * cosz - coszeta * sinz
        prec[1, 1] = -sinzeta * costheta * sinz + coszeta * cosz
        prec[1, 2] = -sinzeta * sintheta
        prec[2, 0] = -sintheta * cosz
        prec[2, 1] = -sintheta * sinz
        prec[2, 2] = costheta

    return prec, psia, wa, ea, xa
