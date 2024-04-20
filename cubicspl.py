def cubicspl(p1, p2, p3, p4):
    """
    Perform cubic splining of an input zero crossing function to find function values.

    Author: David Vallado
    Date: 2 Feb 2004

    Inputs:
        p1, p2, p3, p4: Function values used for splining

    Outputs:
        acu0, acu1, acu2, acu3: Splined polynomial coefficients. acu3*t^3, etc.

    References:
        Vallado, 2013, 559, 1034
    """

    # ------ set up function from C-41 --------
    # det = t1^3*t2^2 + t1^2*t2 + t1*t2^3 - t1^3*t2 - t1^2*t2^3 - t1*t2^2;

    # acu0 = p1;
    # acu1 = ((t2^3-t2^2)*(p2-p1) + (t1^2-t1^3)*(p3-p1) + (t1^3*t2^2-t1^2*t2^3)*(p4-p1)) / det;
    # acu2 = ((t2-t2^3)*(p2-p1) + (t1^3-t1)*(p3-p1) + (t1*t2^3-t1^3*t2)*(p4-p1)) / det;
    # acu3 = ((t2^2-t2)*(p2-p1) + (t1-t1^2)*(p3-p1) + (t1^2*t2-t1*t2^2)*(p4-p1)) / det;

    acu0 = p2
    acu1 = -p1 / 3.0 - 0.5 * p2 + p3 - p4 / 6.0
    acu2 = 0.5 * p1 - p2 + 0.5 * p3
    acu3 = -p1 / 6.0 + 0.5 * p2 - 0.5 * p3 + p4 / 6.0

    return acu0, acu1, acu2, acu3
