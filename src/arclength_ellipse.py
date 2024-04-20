import numpy as np
from elliptic12 import elliptic12

def arclength_ellipse(a, b, theta0=None, theta1=None):
    """
    Calculates the arclength of an ellipse.

    Parameters:
    - a (float): Major axis of the ellipse.
    - b (float): Minor axis of the ellipse.
    - theta0 (float, optional): Starting angle in radians (default: None).
    - theta1 (float, optional): Ending angle in radians (default: None).

    Returns:
    - arclength (float): Length of the arc on the ellipse between theta0 and theta1.
    """

    # Check number of arguments
    if theta0 is None and theta1 is None:
        theta0 = 0
        theta1 = 2 * np.pi
    elif theta0 is None or theta1 is None:
        raise ValueError("Both theta0 and theta1 must be provided if either one is provided.")

    # Default solution for a==b (circles)
    arclength = a * (theta1 - theta0)

    # Ellipses (a<b or a>b)
    if a < b:
        # Theta measured from a axis = semi-MINOR axis
        # Use standard formulation for E(phi,m)
        E1, _ = elliptic12(theta1, 1 - (a / b) ** 2)
        E0, _ = elliptic12(theta0, 1 - (a / b) ** 2)
        arclength = b * (E1 - E0)
    elif a > b:
        # Theta measured from a axis = semi-MAJOR axis
        # Standard formulation will not work ((1-(a/b)^2) < 0); instead use PI/2 - phi and b/a instead of a/b
        E1, _ = elliptic12(np.pi / 2 - theta1, 1 - (b / a) ** 2)
        E0, _ = elliptic12(np.pi / 2 - theta0, 1 - (b / a) ** 2)
        # d(PI/2 - phi)/dphi = -1, so reverse operands in this difference to flip sign:
        arclength = a * (E0 - E1)

    return arclength
