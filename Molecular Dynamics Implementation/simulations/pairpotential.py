from itertools import combinations
from math import sqrt
from typing import Callable

import numpy as np
from numpy import ndarray


def pair_potential(xs: [[float]], potential: Callable[[float, ...], float],
                   potential_args: tuple = ()) -> float:
    """

    Calculates the potential energy of configuration of particles.

    >>> from simulations.potentials import lj_potential as lj

    >>> pair_potential(xs=[[0.0,0.0,0.0]], potential=lj)
    0.0

    >>> pair_potential(xs=[[0.0,0.0,0.0],[0.0,0.0,1.0]], potential=lj)
    0.0

    >>> pair_potential(xs=[[0.0,0.0],[0.0,1.0],[0.0,2.0]],
    ... potential=lj,
    ... potential_args=(1.0, 1.0)) # 2D configuration
    -0.0615234375

    :param xs: positions of the particles
    :param potential: computes the energy of one pair; must be of the form f(x, *args)
    :param potential_args: arguments to pass to the function

    :return: energy of the configuration
    :rtype: float
    """

    energy = 0.0

    for x1, x2 in combinations(xs, 2):  # for statement
        r_squared = 0.0
        for c1, c2 in zip(x1, x2):  # for statement
            r_squared += (c1 - c2) * (c1 - c2)
        r = sqrt(r_squared)  # sqrt imported from math module
        energy += potential(r, *potential_args)

    return energy


def pair_potential_half_vectorised(xs: ndarray, potential: Callable[[float, ...], float],
                                   potential_args: tuple = ()) -> float:
    """

    Calculates the potential energy of configuration of particles.

    >>> from numpy import array
    >>> from simulations.potentials import lj_potential as lj

    >>> pair_potential_half_vectorised(xs=array([[0.0,0.0,0.0]]), potential=lj)
    0.0

    >>> pair_potential_half_vectorised(xs=array([[0.0,0.0,0.0],[0.0,0.0,1.0]]), potential=lj)
    0.0

    >>> pair_potential_half_vectorised(xs=array([[0.0,0.0],[0.0,1.0],[0.0,2.0]]),
    ... potential=lj,
    ... potential_args=(1.0, 1.0)) # 2D configuration
    -0.0615234375

    :param xs: positions of the particles
    :param potential: computes the energy of one pair; must be of the form f(x, *args)
    :param potential_args: arguments to pass to the function

    :return: energy of the configuration
    :rtype: float
    """

    energy = 0.0

    for x1, x2 in combinations(xs, 2):  # for statement
        r = np.linalg.norm(x1-x2)
        energy += potential(r, *potential_args)

    return energy


def pair_potential_vectorised(xs: ndarray, potential: Callable[[ndarray, ...], ndarray],
                              potential_args: tuple = ()) -> float:
    """

    Calculates the potential energy of configuration of particles.

    >>> from numpy import array
    >>> from simulations.potentials import lj_potential_vectorised as lj

    >>> pair_potential_vectorised(xs=array([[0.0,0.0,0.0]]), potential=lj)
    0.0

    >>> pair_potential_vectorised(xs=array([[0.0,0.0,0.0],[0.0,0.0,1.0]]), potential=lj)
    0.0

    >>> pair_potential_vectorised(xs=array([[0.0,0.0],[0.0,1.0],[0.0,2.0]]),
    ... potential=lj,
    ... potential_args=(1.0, 1.0)) # 2D configuration
    -0.0615234375

    :param xs: positions of the particles
    :param potential: computes the energy of one pair; must be of the form f(x, *args)
    :param potential_args: arguments to pass to the function

    :return: energy of the configuration
    :rtype: float
    """

    nparticles, ndim = xs.shape
    # if ndim != 3:
    #     raise ValueError("ndim must be 3")
    left_indices, right_indices = np.triu_indices(nparticles, k=1)
    rij = xs[left_indices] - xs[right_indices]
    dij = np.linalg.norm(rij, axis=1)

    return potential(dij, *potential_args).sum()

def pair_potential_vectorised_spherical(config: ndarray, potential: Callable[[ndarray, ...], ndarray],
                              potential_args: tuple = ()) -> float:
    """
    """
    thetas = config[:, 0]
    phis = config[:, 1]
    num = np.size(thetas)
    xs = np.zeros((num, 3))

    xs[:, 0] = [np.sin(theta) * np.cos(phi) for theta,phi in zip(thetas, phis)]
    xs[:, 1] = [np.sin(theta) * np.sin(phi) for theta,phi in zip(thetas, phis)]
    xs[:, 2] = [np.cos(theta)               for theta     in thetas           ]


    nparticles, ndim = xs.shape
    if ndim != 3:
        raise ValueError("ndim must be 3")
    left_indices, right_indices = np.triu_indices(nparticles, k=1)
    rij = xs[left_indices] - xs[right_indices]
    dij = np.linalg.norm(rij, axis=1)

    return potential(dij, *potential_args).sum()

def pair_potential_pbc(xs: [[float]], potential: Callable[[float, ...], float],
                   potential_args: tuple = ()) -> float:
    """calculate potential energy of configurations using periodical boundary conditon"""
    
    return pair_potential_vectorised(xs, potential, *potential_args)