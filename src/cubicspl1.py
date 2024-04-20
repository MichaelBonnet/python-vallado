from cubic import cubic
from recovpar import recovpar

def cubicspl1(p1, p2, p3, p4):
    """
    Perform cubic splining of an input zero crossing function to find event times.

    Args:
        p1, p2, p3, p4 (float): Function values used for blending.

    Returns:
        minfound (str): Test of success.
        rootf (float): Root for the function.
        funrate (float): Function rate.
    """

    rootf = 0.0
    funrate = 0.0
    minfound = 'n'

    # ---- check for true condition on first two points
    # if (indx == 1) & (sign(p1) ~= sign(p2))
    #     [evt1ctr, evt1, evt2ctr, evt2] = cubicbln ( ev1n, ev2n, timearr, funarr, indx );
    #     event1ctr = event1ctr + evt1ctr;
    #     event1 = [event1;evt1];
    #     event2ctr = event2ctr + evt2ctr;
    #     event2 = [event2;evt2];

    # ------ set up function from C-39 --------
    acu0 = p2
    acu1 = (-p1 + p3) * 0.5
    acu2 = p1 - 2.5 * p2 + 2.0 * p3 - 0.5 * p4
    acu3 = -0.5 * p1 + 1.5 * p2 - 1.5 * p3 + 0.5 * p4

    # --------------- solve roots of this function -------------
    opt = 'U'
    r1r, r1i, r2r, r2i, r3r, r3i = cubic(acu3, acu2, acu1, acu0, opt)

    # ---------- search through roots to locate answers --------
    for indx2 in range(1, 4):
        if indx2 == 1:
            root = r1r
        elif indx2 == 2:
            root = r2r
        elif indx2 == 3:
            root = r3r

        if 0.0 <= root <= 1.0:
            minfound = 'y'
            rootf = root
            # [time] = recovpar(t1,t2,t3,t4, root);
            _, ans = recovpar(p1, p2, p3, p4, root)  # should be 0.0!!!!!!
            # ----- recover the function value derivative
            funrate = 3.0 * acu3 * root ** 2 + 2.0 * acu2 * root + acu1

    # ---- check for true condition on last two points
    # if (sign(p3) ~= sign(p4))
    #     [evt1ctr, evt1, evt2ctr, evt2] = cubicbln ( ev1n, ev2n, timearr, funarr, indx );
    #     event1ctr = event1ctr + evt1ctr;
    #     event1 = [event1;evt1];
    #     event2ctr = event2ctr + evt2ctr;
    #     event2 = [event2;evt2];

    return minfound, rootf, funrate
