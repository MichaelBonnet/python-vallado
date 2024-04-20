import numpy as np
from mag import mag

def angl(vec1, vec2):
    """
    Calculate the angle between two vectors.
    
    Parameters:
        vec1 (array_like): Vector number 1.
        vec2 (array_like): Vector number 2.
        
    Returns:
        float: Angle between the two vectors in the range [-pi, pi]. If the angle
               cannot be determined (due to zero magnitude vectors), returns 999999.1.
    """
    small = 0.00000001
    undefined = 999999.1

    magv1 = mag(vec1)
    magv2 = mag(vec2)

    if magv1 * magv2 > small ** 2:
        temp = np.dot(vec1, vec2) / (magv1 * magv2)
        if abs(temp) > 1.0:
            temp = np.sign(temp) * 1.0
        theta = np.arccos(temp)
    else:
        theta = undefined

    return theta
