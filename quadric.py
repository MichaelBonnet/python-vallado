def quadric(a, b, c, opt):
    """
    Solve for the two roots of a quadratic equation.
    
    This function solves for the two roots of a quadratic equation. There are
    no restrictions on the coefficients, and imaginary results are passed
    out as separate values. The general form is y = ax^2 + bx + c.
    
    Author: David Vallado (719-573-2600), 1 Mar 2001
    Converted to MATLAB by Vallado (719-573-2600), 3 Dec 2002

    Inputs:
        a   - coefficient of x squared term
        b   - coefficient of x term
        c   - constant
        opt - option for output
              'I' - all roots including imaginary
              'R' - only real roots
              'U' - only unique real roots (no repeated)
    
    Outputs:
        r1r - real portion of root 1
        r1i - imaginary portion of root 1
        r2r - real portion of root 2
        r2i - imaginary portion of root 2
    
    """
    # Implementation
    small = 0.00000001
    r1r, r1i, r2r, r2i = 0.0, 0.0, 0.0, 0.0
    
    discrim = b * b - 4.0 * a * c
    
    # Real roots
    if abs(discrim) < small:
        r1r = -b / (2.0 * a)
        r2r = r1r
    else:
        if abs(a) < small:
            r1r = -c / b
        else:
            if discrim > 0.0:
                r1r = (-b + discrim ** 0.5) / (2.0 * a)
                r2r = (-b - discrim ** 0.5) / (2.0 * a)
            else:
                # Complex roots
                if opt == 'I':
                    r1r = -b / (2.0 * a)
                    r2r = r1r
                    r1i = (abs(-discrim) ** 0.5) / (2.0 * a)
                    r2i = -(abs(-discrim) ** 0.5) / (2.0 * a)
                else:
                    r1r = 99999.9
                    r2r = 99999.9
    
    return r1r, r1i, r2r, r2i
