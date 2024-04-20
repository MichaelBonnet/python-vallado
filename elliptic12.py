import numpy as np

def elliptic12(u, m, tol=np.finfo(float).eps):
    """
    Evaluate the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function.

    Parameters:
    u (array_like): Phase in radians.
    m (array_like): Module (0 < m < 1).
    tol (float, optional): Tolerance. Default is machine epsilon.

    Returns:
    F (ndarray): Incomplete Elliptic Integral of the First Kind.
    E (ndarray): Incomplete Elliptic Integral of the Second Kind.
    Z (ndarray): Jacobi Zeta Function.

    Uses the method of the Arithmetic-Geometric Mean and Descending Landen Transformation.
    """
    if np.ndim(u) == 0:
        u = np.array([u])
    if np.ndim(m) == 0:
        m = np.array([m])
    
    if len(u) != len(m):
        raise ValueError("U and M must be the same size.")
    
    F = np.zeros_like(u)
    E = np.zeros_like(u)
    Z = np.zeros_like(u)

    if np.any(np.logical_or(m < 0, m > 1)):
        raise ValueError("M must be in the range 0 <= M <= 1.")
    
    m = np.asarray(m)
    u = np.asarray(u)

    # cdav change for small eccentricities
    m[m < 1e-7] = 1e-7

    I = np.where(np.logical_and(m != 1, m != 0))[0]
    if len(I) > 0:
        mu, J, K = np.unique(m[I], return_index=True, return_inverse=True)
        mumax = len(mu)
        signU = np.sign(u[I])

        chunk = 7
        a = np.zeros((chunk, mumax))
        c = np.copy(a)
        b = np.copy(a)
        a[0, :] = 1
        c[0, :] = np.sqrt(mu)
        b[0, :] = np.sqrt(1 - mu)
        n = np.zeros(mumax, dtype=np.uint32)
        i = 0
        while np.any(np.abs(c[i, :]) > tol):
            i += 1
            if i >= a.shape[0]:
                a = np.concatenate((a, np.zeros((2, mumax))))
                b = np.concatenate((b, np.zeros((2, mumax))))
                c = np.concatenate((c, np.zeros((2, mumax))))
            a[i, :] = 0.5 * (a[i-1, :] + b[i-1, :])
            b[i, :] = np.sqrt(a[i-1, :] * b[i-1, :])
            c[i, :] = 0.5 * (a[i-1, :] - b[i-1, :])
            in_indices = np.where(np.logical_and(np.abs(c[i, :]) <= tol, np.abs(c[i-1, :]) > tol))[0]
            if len(in_indices) > 0:
                n[in_indices] = i - 1
        
        mmax = len(I)
        mn = np.max(n)
        phin = np.zeros(mmax)
        C = np.zeros(mmax)
        Cp = np.zeros(mmax)
        e = np.zeros(mmax, dtype=np.uint32)
        phin[:] = signU * u[I]
        i = 0
        while i < mn:
            i += 1
            in_indices = np.where(n[K] > i)[0]
            if len(in_indices) > 0:
                phin[in_indices] = np.arctan(b[i, K[in_indices]] / a[i, K[in_indices]] * np.tan(phin[in_indices])) \
                                    + np.pi * (np.ceil(phin[in_indices] / np.pi - 0.5)) + phin[in_indices]
                e[in_indices] = 2 ** (i - 1)
                C[in_indices] += e[in_indices[0]] * c[i, K[in_indices]]
                Cp[in_indices] += c[i + 1, K[in_indices]] * np.sin(phin[in_indices])
        
        Ff = phin / (a[mn, K] * e * 2)
        F[I] = Ff * signU
        Z[I] = Cp * signU
        E[I] = (Cp + (1 - 1/2 * C) * Ff) * signU

    # Special cases: m == {0, 1}
    m0 = np.where(m == 0)[0]
    if len(m0) > 0:
        F[m0] = u[m0]
        E[m0] = u[m0]

    m1 = np.where(m == 1)[0]
    um1 = np.abs(u[m1])
    if len(m1) > 0:
        N = np.floor((um1 + np.pi / 2) / np.pi)
        M = np.where(um1 < np.pi / 2)[0]
        
        F[m1[M]] = np.log(np.tan(np.pi / 4 + u[m1[M]] / 2))
        F[m1[um1 >= np.pi / 2]] = np.inf * np.sign(u[m1[um1 >= np.pi / 2]])

        E[m1] = ((-1) ** N * np.sin(um1) + 2 * N) * np.sign(u[m1])

        Z[m1] = (-1) ** N * np.sin(u[m1])

    return F, E, Z
