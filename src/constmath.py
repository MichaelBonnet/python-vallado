"""
Constants for mathematical operations.

Author:
    David Vallado
    Phone: 719-573-2600
    Date: 2 Apr 2007

Locals:
    small: A small value
    infinite: Representation of infinity
    undefined: Representation of undefined value
"""

import math

small = 1.0e-10

infinite = 999999.9
undefined = 999999.1

# -------------------------  mathematical  --------------------
rad = 180.0 / math.pi # Radian conversion factor
twopi = 2.0 * math.pi # 2 * pi
halfpi = math.pi * 0.5 # 0.5 * pi

# -------------------------  conversions  ---------------------
ft2m = 0.3048 # Feet to meters conversion factor
mile2m = 1609.344 # Miles to meters conversion factor
nm2m = 1852 # Nautical miles to meters conversion factor
mile2ft = 5280 # Miles to feet conversion factor
mileph2kmph = 0.44704 # Miles per hour to kilometers per hour conversion factor
nmph2kmph = 0.5144444 # Nautical miles per hour to kilometers per hour conversion factor
