def azl2radc(az, el, lat, lst):
    """
    azl2radc

    This function finds the rtasc decl values given the az-el.

    Args:
    - az (float): Azimuth angle.
    - el (float): Elevation angle.
    - lat (float): Latitude.
    - lst (float): Local Sidereal Time.

    Returns:
    - rtasc (float): Right Ascension.
    - decl (float): Declination.
    """

    import math

    rad = 180.0 / math.pi

    decl = math.asin(math.sin(el) * math.sin(lat) + math.cos(el) * math.cos(lat) * math.cos(az))

    slha1 = -(math.sin(az) * math.cos(el) * math.cos(lat)) / (math.cos(decl) * math.cos(lat))
    clha1 = (math.sin(el) - math.sin(lat) * math.sin(decl)) / (math.cos(decl) * math.cos(lat))

    lha1 = math.atan2(slha1, clha1)
    # print(' lha1 %13.7f \n' % (lha1 * rad))

    # alt approach
    slha2 = -(math.sin(az) * math.cos(el)) / (math.cos(decl))
    clha2 = (math.cos(lat) * math.sin(el) - math.sin(lat) * math.cos(el) * math.cos(az)) / (math.cos(decl))

    lha2 = math.atan2(slha2, clha2)
    # print(' lha2 %13.7f \n' % (lha2 * rad))

    rtasc = lst - lha1

    return rtasc, decl
