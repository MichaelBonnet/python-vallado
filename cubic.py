from quadric import quadric

def cubic(a3, b2, c1, d0, opt):
    """
    Solve for the three roots of a cubic equation.

    Parameters:
    a3 : float
        Coefficient of x cubed term.
    b2 : float
        Coefficient of x squared term.
    c1 : float
        Coefficient of x term.
    d0 : float
        Constant term.
    opt : str
        Option for output:
        'I' for all roots including imaginary.
        'R' for only real roots.
        'U' for only unique real roots (no repeated).

    Returns:
    r1r : float
        Real portion of root 1.
    r1i : float
        Imaginary portion of root 1.
    r2r : float
        Real portion of root 2.
    r2i : float
        Imaginary portion of root 2.
    r3r : float
        Real portion of root 3.
    r3i : float
        Imaginary portion of root 3.
    """
    # --------------------  implementation   ----------------------
    import math

    rad = 180.0 / math.pi
    onethird = 1.0 / 3.0
    small = 0.0000000001
    r1r = 0.0
    r1i = 0.0
    r2r = 0.0
    r2i = 0.0
    r3r = 0.0
    r3i = 0.0

    if abs(a3) > small:
        # ----------- force coefficients into std form ----------------
        p = b2 / a3
        q = c1 / a3
        r = d0 / a3

        a3 = onethird * (3.0 * q - p * p)
        b2 = (1.0 / 27.0) * (2.0 * p * p * p - 9.0 * p * q + 27.0 * r)

        delta = (a3 * a3 * a3 / 27.0) + (b2 * b2 * 0.25)

        # ------------------ use cardans formula ----------------------
        if delta > small:
            temp1 = (-b2 * 0.5) + math.sqrt(delta)
            temp2 = (-b2 * 0.5) - math.sqrt(delta)
            temp1 = math.copysign(abs(temp1) ** onethird, temp1)
            temp2 = math.copysign(abs(temp2) ** onethird, temp2)
            r1r = temp1 + temp2 - p * onethird

            if opt == 'I':
                r2r = -0.5 * (temp1 + temp2) - p * onethird
                r2i = -0.5 * math.sqrt(3.0) * (temp1 - temp2)
                r3r = -0.5 * (temp1 + temp2) - p * onethird
                r3i = -r2i
            else:
                r2r = 99999.9
                r3r = 99999.9
        else:
            # --------------- evaluate zero point ---------------------
            if abs(delta) < small:
                r1r = -2.0 * math.copysign(abs(b2 * 0.5) ** onethird, b2) - p * onethird
                r2r = math.copysign(abs(b2 * 0.5) ** onethird, b2) - p * onethird
                r3r = r2r
            else:
                # ------------ use trigonometric identities -----------
                e0 = 2.0 * math.sqrt(-a3 * onethird)
                cosphi = (-b2 / (2.0 * math.sqrt(-a3 * a3 * a3 / 27.0)))
                sinphi = math.sqrt(1.0 - cosphi * cosphi)
                phi = math.atan2(sinphi, cosphi)
                if phi < 0.0:
                    phi = phi + 2.0 * math.pi
                r1r = e0 * math.cos(phi * onethird) - p * onethird
                r2r = e0 * math.cos(phi * onethird + 120.0 / rad) - p * onethird
                r3r = e0 * math.cos(phi * onethird + 240.0 / rad) - p * onethird
    else:
        r1r, r1i, r2r, r2i = quadric(b2, c1, d0, opt)
        r3r = 99999.9
        r3i = 99999.9

    return r1r, r1i, r2r, r2i, r3r, r3i