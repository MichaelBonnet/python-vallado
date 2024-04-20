"""
Constants for various astrodynamic operations.

Author: David Vallado (719-573-2600), 2 Apr 2007
"""

import numpy as np

# EGM-08 constants used here
re = 6378.1363  # Earth's equatorial radius (km)
flat = 1.0 / 298.257223563 # Earth's flattening
earthrot = 7.292115e-5  # Earth's rotation rate (rad/s) old 7.29211514670698e-05
mu = 398600.4415  # Earth's gravitational parameter (km^3/s^2)
mum = 3.986004415e14  # Earth's gravitational parameter (m^3/s^2)

# derived constants from the base values
eccearth = np.sqrt(2.0 * flat - flat ** 2) # Earth's eccentricity
eccearthsqrd = eccearth ** 2 # Square of Earth's eccentricity

nm2m = 1852.0 # nautical miles to meters
ft2m = 0.3048 # feet to meters

renm = re / nm2m # Earth's equatorial radius in nautical miles
reft = re * 1000.0 / ft2m # Earth's equatorial radius in feet

tusec = np.sqrt(re ** 3 / mu) # Time unit in seconds
tumin = tusec / 60.0 # Time unit in minutes
tuday = tusec / 86400.0 # Time unit in days
tudaysid = tusec / 86164.090524 # Time unit in sidereal days

omegaearthradptu = earthrot * tusec # Earth's rotation rate in radian per time unit
omegaearthradpmin = earthrot * 60.0 # Earth's rotation rate in radian per minute

velkmps = np.sqrt(mu / re) # Earth's circular velocity (km/s)
velftps = velkmps * 1000.0 / ft2m # Earth's circular velocity (ft/s)
velradpmin = velkmps * 60.0 / re # Earth's circular velocity in radian per minute
# for afspc
# velkmps1 = velradpmin*6378.135/60.0   7.90537051051763
# mu1 = velkmps*velkmps*6378.135        3.986003602567418e+005
degpsec = (180.0 / np.pi) / tusec # Degrees per second
radpday = 2.0 * np.pi * 1.002737909350795 # Radians per day

speedoflight = 299792.458  # Speed of light (km/s)
au = 149597870.7  # Astronomical unit (km)
earth2moon = 384400.0  # Distance from Earth to Moon (km)
moonradius = 1738.0  # Moon's radius (km)
sunradius = 696000.0  # Sun's radius (km)

masssun = 1.9891e30 # Mass of the Sun (kg)
massearth = 5.9742e24 # Mass of the Earth (kg)
massmoon = 7.3483e22 # Mass of the Moon (kg)

musun = 1.32712428e11 # Gravitational parameter of the Sun (km^3/s^2)
mumoon = 4902.799 # Gravitational parameter of the Moon (km^3/s^2)
