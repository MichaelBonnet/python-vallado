import numpy as np
from mag import mag
from angl import angl
from rv2radec import rv2radec
from rv2coe import rv2coe
from matmult import matmult, matvecmult
from doubler import doubler

# Function to solve the problem of orbit determination using three optical sightings
def anglesdr(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3,
             rsite1, rsite2, rsite3, re, mu):
    # Constants
    rad = 180.0 / np.pi
    magr1in = 2.0 * re
    magr2in = 2.04 * re
    direct = 'y'
    tol = 1e-8 * re  # km
    pctchg = 0.005

    # Convert days to seconds
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0
    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0

    # Form line of sight vectors
    los1 = np.array([np.cos(decl1) * np.cos(rtasc1), np.cos(decl1) * np.sin(rtasc1), np.sin(decl1)])
    los2 = np.array([np.cos(decl2) * np.cos(rtasc2), np.cos(decl2) * np.sin(rtasc2), np.sin(decl2)])
    los3 = np.array([np.cos(decl3) * np.cos(rtasc3), np.cos(decl3) * np.sin(rtasc3), np.sin(decl3)])

    # Other initializations
    magr1old = 99999.9
    magr2old = 99999.9
    magrsite1 = mag(rsite1)
    magrsite2 = mag(rsite2)
    magrsite3 = mag(rsite3)
    cc1 = 2.0 * np.dot(los1, rsite1)
    cc2 = 2.0 * np.dot(los2, rsite2)
    ktr = 0

    # Main loop for double-r algorithm
    while (np.abs(magr1in - magr1old) > tol or np.abs(magr2in - magr2old) > tol):
        ktr += 1
        print('%2i ' % ktr)
        # Call doubler function
        [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                                                   los1, los2, los3, rsite1, rsite2, rsite3,
                                                                   tau12, tau32, direct, re, mu)
        # Check intermediate status
        f = 1.0 - a / magr2 * (1.0 - np.cos(deltae32))
        g = tau32 - np.sqrt(a ** 3 / mu) * (deltae32 - np.sin(deltae32))
        v2 = (r3 - f * r2) / g
        [p, a, ecc, incl, omega, argp, nu, m] = rv2coe(r2, v2, re, mu)
        print('coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f' % (
            p, a, ecc, incl * rad, omega * rad, argp * rad, nu * rad, m * rad))

        # Recalculate f1 and f2 with r1 = r1 + delta r1
        magr1o = magr1in
        deltar1 = pctchg * magr1in
        magr1in = magr1in + deltar1
        [r2, r3, f1delr1, f2delr1, q2, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in,
                                                                              magr2in, los1, los2, los3, rsite1,
                                                                              rsite2, rsite3, tau12, tau32, direct,
                                                                              re, mu)
        pf1pr1 = (f1delr1 - f1) / deltar1
        pf2pr1 = (f2delr1 - f2) / deltar1

        # Recalculate f1 and f2 with r2 = r2 + delta r2
        magr1in = magr1o
        magr2o = magr2in
        deltar2 = pctchg * magr2in
        magr2in = magr2in + deltar2
        [r2, r3, f1delr2, f2delr2, q3, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in,
                                                                              magr2in, los1, los2, los3, rsite1,
                                                                              rsite2, rsite3, tau12, tau32, direct,
                                                                              re, mu)
        pf1pr2 = (f1delr2 - f1) / deltar2
        pf2pr2 = (f2delr2 - f2) / deltar2

        # Calculate updates
        delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2
        delta1 = pf2pr2 * f1 - pf1pr2 * f2
        delta2 = pf1pr1 * f2 - pf2pr1 * f1
        deltar1 = -delta1 / delta
        deltar2 = -delta2 / delta
        magr1old = magr1in
        magr2old = magr2in

        # Limit the amount of correction
        if np.abs(deltar1) > magr1in * pctchg:
            print('%11.7f ' % deltar1)
        if np.abs(deltar2) > magr2in * pctchg:
            print('%11.7f ' % deltar2)

        magr1in = magr1in + deltar1
        magr2in = magr2in + deltar2

        print('qs %11.7f  %11.7f  %11.7f' % (q1, q2, q3))
        print('magr1o %11.7f delr1 %11.7f magr1 %11.7f %11.7f' % (magr1o, deltar1, magr1in, magr1old))
        print('magr2o %11.7f delr2 %11.7f magr2 %11.7f %11.7f' % (magr2o, deltar2, magr2in, magr2old))
        print('=============================================== \n')

    # Needed to set r2 properly since the last one was moving r2
    [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                                               los1, los2, los3, rsite1, rsite2, rsite3,
                                                               tau12, tau32, direct, re, mu)

    f = 1.0 - a / magr2 * (1.0 - np.cos(deltae32))
    g = tau32 - np.sqrt(a ** 3 / mu) * (deltae32 - np.sin(deltae32))
    v2 = (r3 - f * r2) / g

    return r2, v2
