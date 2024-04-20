from cubicspl import cubicspl
from cubic import cubic

def ap2kp(apin):
    """
    Convert ap to kp using cubic splines.

    Author: David Vallado, 719-573-2600, 4 Aug 2005.

    Inputs:
    - apin: ap

    Outputs:
    - kpout: kp

    References:
    - Vallado, 2004, 899-901

    kpout = ap2kp(apin)
    """
    ap = [-0.00001, -0.001, 0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179,
          207, 236, 300, 400, 900]
    kp = [-0.66666667, -0.33333, 0.0, 0.33333, 0.66667, 1.0, 1.33333, 1.66667, 2.0, 2.33333, 2.66667,
          3.0, 3.33333, 3.66667, 4.0, 4.33333, 4.66667, 5.0, 5.33333, 5.66667,
          6.0, 6.33333, 6.66667, 7.0, 7.33333, 7.66667, 8.0, 8.33333, 8.66667, 9.0, 9.33333]

    kpout = 0.0
    t1 = 0.0
    t2 = 0.0

    i = 0
    while i < 30 and apin > ap[i]:
        i += 1

    if i > 3:
        p1 = kp[i - 2]
        p2 = kp[i - 1]
        p3 = kp[i]
        p4 = kp[i + 1]
        ac0, ac1, ac2, ac3 = cubicspl(p1, p2, p3, p4, t1, t2)

        p1 = ap[i - 2]
        p2 = ap[i - 1]
        p3 = ap[i]
        p4 = ap[i + 1]
        kac0, kac1, kac2, kac3 = cubicspl(p1, p2, p3, p4, t1, t2)

        r1r, _, r2r, _, r3r, _ = cubic(kac3, kac2, kac1, kac0 - apin, 'R')

        r = [r1r, r2r, r3r]

        if 0.00 <= r[0] <= 1.001:
            apval = r[0]
        elif 0.00 <= r[1] <= 1.001:
            apval = r[1]
        elif 0.00 <= r[2] <= 1.001:
            apval = r[2]
        else:
            apval = 0.0
            print('error in root %11.7f %3i %11.7f %11.7f %11.7f' % (apin, i, r[0], r[1], r[2]))

        kpout = ac3 * apval**3 + ac2 * apval**2 + ac1 * apval + ac0

    return kpout