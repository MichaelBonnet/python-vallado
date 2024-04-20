import numpy as np

# Function to perform matrix multiplication
def matmult(a, b):
    return np.dot(a, b)

def matvecmult(mat, vec, size_of):
    """
    ----------------------------------------------------------------------------
    
                                   function matvecmult
    
      this function multiplies a matrix by a vector.
    
      author        : david vallado                  719-573-2600    4 jun 2002
    
      revisions
                    -
    
      inputs          description                    range / units
        mat         - matrix (square)
        vec         - vector
        size        - dimension of matrix
    
      outputs       :
        outvec      - unit vector
    
      locals        :
        i,j         - index
    
      coupling      :
        none
    
      [outvec] = matvecmult ( mat, vec, sizeof );
    ----------------------------------------------------------------------------
    """
    # -------------------------  implementation   -----------------
    outvec = [0.0] * size_of
    for i in range(size_of):
        for j in range(size_of):
            outvec[i] += mat[i][j] * vec[j]
    return outvec
