import numpy as np
from numpy import ndarray


def lj_potential(r: float, epsilon: float = 1.0, sigma: float = 1.0) -> float:
    """
    compute the Lennard Jones potential at particle separation r,

    V_LJ = 4 epsilon ( (sigma/r)^12 - (sigma/r)^6 )

    >>> lj_potential(2**(1/6), 1,1)
    -1.0

    >>> lj_potential(-1)
    Traceback (most recent call last):
    ValueError: r should be greater than 0, not -1

    >>> lj_potential('1')
    Traceback (most recent call last):
    TypeError: r should be float, not str

    """

    if not isinstance(r, (int, float)):
        raise TypeError(f'r should be float, not {type(r).__name__}')
    elif r <= 0:
        raise ValueError(f'r should be greater than 0, not {r}')

    r6 = (sigma / r) ** 2
    r6 *= r6 * r6
    return 4 * epsilon * r6 * (r6 - 1)


def lj_potential_vectorised(rs: ndarray, epsilon: float = 1.0, sigma: float = 1.0) -> ndarray:
    """
    compute the Lennard Jones potential at particle separation r,

    V_LJ = 4 epsilon ( (sigma/r)^12 - (sigma/r)^6 )

    >>> from numpy import array
    >>> lj_potential_vectorised(array([1.0, 2**(1/6), 2.0]), epsilon=1.0, sigma=1.0)  #doctest: +ELLIPSIS
    array([ 0.        , -1.        , -0.0615...])

    >>> lj_potential_vectorised(-1)
    Traceback (most recent call last):
    TypeError: rs should be ndarray, not int

    >>> lj_potential_vectorised(array([1.0, 2**(1/6), -1.0]))
    Traceback (most recent call last):
    ValueError: all r must be positive, but min(r) = -1.0

    >>> lj_potential_vectorised(array([1.0, 2**(1/6), 2.0]), epsilon=-1.0, sigma=1.0)
    Traceback (most recent call last):
    ValueError: both epsilon and sigma must be positive, not (-1.0, 1.0)
    """

    if not isinstance(rs, ndarray):
        raise TypeError(f'rs should be ndarray, not {type(rs).__name__}')
    if np.any(rs <= 0):
        raise ValueError(f'all r must be positive, but min(r) = {rs.min()}')
    if epsilon <= 0 or sigma <= 0:
        raise ValueError(f'both epsilon and sigma must be positive, not ({epsilon}, {sigma})')

    r6 = (sigma / rs) ** 2
    r6 *= r6 * r6
    return 4 * epsilon * r6 * (r6 - 1)

def coloumb_potential_vectorised(rs: ndarray) -> ndarray:
    """
    compute the coloumb potential at particle separation r,

    V_CB = 1 / rs^2

    """

    if not isinstance(rs, ndarray):
        raise TypeError(f'rs should be ndarray, not {type(rs).__name__}')
    if np.any(rs <= 0):
        raise ValueError(f'all r must be positive, but min(r) = {rs.min()}')

    return 1 / rs**2


def cut_off_potential(rs: ndarray) -> ndarray:
    """cut off potential at particle separation r,
    
    """
    rs = np.array(rs)
    rc = 1.5 # cutoff r
    beta = rc / 1
    alpha = (27 / 4) * (beta ** 2 / (1 - beta ** 2) ** 3)
    res = []
    for r in rs:
        if r <= rc:
            res.append( - alpha * (1 / r **2 - 1) * (beta ** 2 / r ** 2 - 1) ** 2)
        else:
            res.append(0)

    return np.array(res)