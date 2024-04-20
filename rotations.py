import numpy as np

def rot1(vec, xval):
    """
    Perform a rotation about the 1st axis.

    Parameters:
    - vec : list or tuple
        Input vector.
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outvec : list
        Vector result.
    """
    temp = vec[2]
    c = np.cos(xval)
    s = np.sin(xval)

    outvec = [0, 0, 0]
    outvec[2] = c * vec[2] - s * vec[1]
    outvec[1] = c * vec[1] + s * temp
    outvec[0] = vec[0]

    return outvec


def rot1mat(xval):
    """
    Set up a rotation matrix for an input angle about the first axis.

    Parameters:
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outmat : list of lists
        Matrix result.
    """
    c = np.cos(xval)
    s = np.sin(xval)

    outmat = [[1.0, 0.0, 0.0],
              [0.0, c, s],
              [0.0, -s, c]]

    return outmat


def rot2(vec, xval):
    """
    Perform a rotation about the 2nd axis.

    Parameters:
    - vec : list or tuple
        Input vector.
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outvec : list
        Vector result.
    """
    temp = vec[2]
    c = np.cos(xval)
    s = np.sin(xval)

    outvec = [0, 0, 0]
    outvec[2] = c * vec[2] + s * vec[0]
    outvec[0] = c * vec[0] - s * temp
    outvec[1] = vec[1]

    return outvec


def rot2mat(xval):
    """
    Set up a rotation matrix for an input angle about the second axis.

    Parameters:
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outmat : list of lists
        Matrix result.
    """
    c = np.cos(xval)
    s = np.sin(xval)

    outmat = [[c, 0.0, -s],
              [0.0, 1.0, 0.0],
              [s, 0.0, c]]

    return outmat


def rot3(vec, xval):
    """
    Perform a rotation about the 3rd axis.

    Parameters:
    - vec : list or tuple
        Input vector.
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outvec : list
        Vector result.
    """
    temp = vec[1]
    c = np.cos(xval)
    s = np.sin(xval)

    outvec = [0, 0, 0]
    outvec[1] = c * vec[1] - s * vec[0]
    outvec[0] = c * vec[0] + s * temp
    outvec[2] = vec[2]

    return outvec


def rot3mat(xval):
    """
    Set up a rotation matrix for an input angle about the third axis.

    Parameters:
    - xval : float
        Angle of rotation in radians.

    Returns:
    - outmat : list of lists
        Matrix result.
    """
    c = np.cos(xval)
    s = np.sin(xval)

    outmat = [[c, s, 0.0],
              [-s, c, 0.0],
              [0.0, 0.0, 1.0]]

    return outmat
