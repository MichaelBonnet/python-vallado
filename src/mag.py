import numpy as np

def mag(vec):
    """
    Calculate the magnitude of a vector.
    
    Parameters:
        vec (array_like): Input vector.
        
    Returns:
        float: Magnitude of the vector.
    """
    temp = vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2

    if abs(temp) >= 1.0e-16:
        mag = np.sqrt(temp)
    else:
        mag = 0.0

    return mag
