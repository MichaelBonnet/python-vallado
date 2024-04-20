def polarm(xp, yp, ttt, opt):
    """
    Calculates the transformation matrix that accounts for polar motion.
    
    Both the 1980 and 2000 theories are handled. Note that the rotation 
    order is different between 1980 and 2000.
    
    Author: David Vallado (719-573-2600), 25 Jun 2002
    
    Inputs:
        xp : float
            Polar motion coefficient in radians
        yp : float
            Polar motion coefficient in radians
        ttt : float
            Julian centuries of TT (00 theory only)
        opt : str
            Method option: '01', '02', '80'
    
    Outputs:
        pm : list of lists
            Transformation matrix for ECEF - PEF
    
    References:
        Vallado 2004, 207-209, 211, 223-224
    """
    import math
    
    cosxp = math.cos(xp)
    sinxp = math.sin(xp)
    cosyp = math.cos(yp)
    sinyp = math.sin(yp)

    pm = [[0.0] * 3 for _ in range(3)]

    if opt == "80":
        pm[0][0] = cosxp
        pm[0][1] = 0.0
        pm[0][2] = -sinxp
        pm[1][0] = sinxp * sinyp
        pm[1][1] = cosyp
        pm[1][2] = cosxp * sinyp
        pm[2][0] = sinxp * cosyp
        pm[2][1] = -sinyp
        pm[2][2] = cosxp * cosyp

        # a1 = rot2mat(xp);
        # a2 = rot1mat(yp);
        # pm = a2*a1;          
        # Approximate matrix using small angle approximations
        # pm[0][0] =  1.0
        # pm[1][0] =  0.0
        # pm[2][0] =  xp
        # pm[0][1] =  0.0
        # pm[1][1] =  1.0
        # pm[2][1] = -yp
        # pm[0][2] = -xp
        # pm[1][2] =  yp
        # pm[2][2] =  1.0
    else:
        convrt = math.pi / (3600.0 * 180.0)
        # approximate sp value in rad
        sp = -47.0e-6 * ttt * convrt
        print('xp {:16.14f}, {:16.14f} sp {:16.14g}'.format(xp, yp, sp))
        cossp = math.cos(sp)
        sinsp = math.sin(sp)

        # print(' sp  {:14.11f} mas'.format(sp/convrt))

        # form the matrix
        pm[0][0] = cosxp * cossp
        pm[0][1] = -cosyp * sinsp + sinyp * sinxp * cossp
        pm[0][2] = -sinyp * sinsp - cosyp * sinxp * cossp
        pm[1][0] = cosxp * sinsp
        pm[1][1] = cosyp * cossp + sinyp * sinxp * sinsp
        pm[1][2] = sinyp * cossp - cosyp * sinxp * sinsp
        pm[2][0] = sinxp
        pm[2][1] = -sinyp * cosxp
        pm[2][2] = cosyp * cosxp

        # a1 = rot1mat(yp);
        # a2 = rot2mat(xp);
        # a3 = rot3mat(-sp);
        # pm = a3*a2*a1;

    return pm
