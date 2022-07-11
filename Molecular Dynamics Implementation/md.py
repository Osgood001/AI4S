import numpy as np
from typing import List, Callable
from itertools import product, combinations
import matplotlib.pyplot as plt
import json
from numpy import ndarray
from simulations import pair_potential, potentials, pair_potential_pbc, cut_off_potential
from tqdm import tqdm

r = 0.5

# generate ordered configuration
def ordered_config(N: int, r: float, L: float):
    """Generate ordered configuration in a L*L size box with N particles of radius r
    return coordinates, number_placed and cell_list
    """

    ...

# generate random configuration
def random_config(N: int, r: float, L: float):
    """generate random configuration in a L*L size box with N particles of radius r
    note that L must be integer multiple of diameter (2 * r)
    # >>> L = 5
    # >>> r = 0.5
    # >>> N = 10
    # >>> coord, num_placed, cell_list = random_config(N, r, L)
    # >>> len(coord)
    # 10

    """
    size = int(L / (2 * r))  # don't add one

    grid_index = list(product(range(size), range(size)))# split a box into grid
    num_placed = 0                                      # number of particles been placed in the grid
    res = {}
    while num_placed < N and grid_index:
        rand_index = np.random.randint(len(grid_index)) # choose random index in the grid
        i,j = grid_index[rand_index]
        del grid_index[rand_index]                      # remove the index from the grid index

        temp = 0
        attemp_iter = 6 * 1e2
        while temp < attemp_iter: 
            left, right = 2 * r * j, 2 * r * (j + 1)
            down, up = 2 * r * i , 2 * r * (i + 1)

            # create a random coorddinate in the box
            x = np.random.rand() * (right - left) + left
            y = np.random.rand() * (up - down) + down

            nearby_index = nearby(str([i, j]), size)
            valid = 1 # if the choosen random coordinate is valid
            for item in nearby_index:

                if valid == 0: break
                try:
                    coordinates = res[str(item)]
                    for x0, y0 in coordinates:
                        if (x - x0)**2 + (y - y0)**2 <= (2 * r)**2:
                            valid = 0
                            break
                except KeyError:
                    continue
            if valid: # vaid for all particles nearby
                item = [i,j]
                num_placed += 1
                # timing
                timing(num_placed, N)
                try:
                    res[str(item)]
                    res[str(item)].append([x, y])
                except KeyError:
                    res[str(item)] = [[x, y]]
                break
            temp += 1
    
    answer = []
    for key, value in res.items():
        for item in value:
            answer.append(item)
    return answer, num_placed, res

# run monte carol simulation
def monte_carol(cell: dict, T: float, u0: float, m: int=100):
    """return a sequence of energy  
    [u0, u1, u2, ... , um] for MC simulation
    
    cell: the dict we obtained in the generation period, 
    this data structure will be used generally as the configuration

    T: the temperature of the monte carol simulation

    u0: the energy of the configuration recorded by cell
     (can be obtained by previous code)

    m: step of valid MC needed, 100 as default
    """
    visualize(cell)
    n = 0 # valid step of MC simulation
    energy = [u0] # energy sequence
    while n < m:
        # timing
        timing(n, m)

        REJECT = 1 # the status variable of whether the pertubation is rejected
        while REJECT:
            cell_key = list(cell.keys())                            # the keys list of the cell
            cell_index = cell_key[np.random.randint(len(cell))]     # choosen cell
            nearby_index = nearby(cell_index)                       # nearby cells
            particle_index = np.random.randint(len(cell[cell_index]))
            center = cell[cell_index][particle_index]               # chosen particle as center
            pert_center = pertubation(center)                       # pertubation of the chosen particle

            # find the pertubated box index and the pertubated surroundings
            pert_cell_index = [int(pert_center[0] / (2 * r)), int(pert_center[1] / (2 * r))]
            pert_nearby_index = nearby(str(pert_cell_index))

            u1 = center_potential(center, nearby_index, cell)
            u2 = center_potential(pert_center, pert_nearby_index, cell)
        
            if metropolis(u2 - u1, T):
                del cell[cell_index][particle_index]
                if len(cell[cell_index]) < 1: del cell[cell_index]
                try: 
                    cell[str(pert_cell_index)].append(pert_center)  # update the particle position
                except KeyError:
                    cell[str(pert_cell_index)] = [pert_center]  # update the particle position
                energy.append(energy[-1] + u2 - u1)
                REJECT = 0
                n += 1
    visualize(cell)
    return energy

# molecular dynamics: velocity-verlet algorithm
def velocity_verlet(coord: np.ndarray, velocity: np.ndarray, potential: Callable, period: float, step: int, box_length: float):
    """Velocity-Verlet algorithm 
    
    Args:
        coord: coordinates of all particles
        velocity: velocity of particles
        potential: pairwise potential (refer to lj potential implemented before)
        period: the length of one time step
        step: the steps needed for the algorithm to stop

    Returns:
        final coordinates
        final velocity
        energy sequence
    """
    L = box_length # length of the box
    t = period # time interval
    mass = 1 # set default mass to be 1
    cR = np.array(coord) # current coordinates R
    cV = velocity # current velocity
    cA = force(cR) / mass # current acReleration
    nR = update_coord(cR, cV, cA, t, L) # new coordinates
    nA = force(nR) / mass # new acceleration
    nV = update_velocity(cV, cA, nA, t) # new velocity
    energy = [pair_potential_pbc(cR, cut_off_potential)] # energy sequence ! pbc needed

    for i in tqdm(range(step)):
        cR = nR
        cV = nV
        cA = nA

        nR = update_coord(cR, cV, cA, t, L)
        nA = force(nR) / mass
        nV = update_velocity(cV, cA, nA, t)

        energy.append(pair_potential_pbc(cR, cut_off_potential))

    return cR, cV, energy


# calculate the force act on all particles
def force(cR: np.ndarray, epsilon=1, alpha=1,rc=1.25,sigma=1) -> float:
    """ calculatet the force of the cut off pootential
    
    analytic discription is Fr = epsilon * alpha * [(4 * rc**2 - 2 * sigma**2) * r**-3 \
                                             - (4 * rc**4 + 2 * sigma**2 * rc**2) * r**-5 
                                             + 4 * sigma**2 * rc**4 * r**-7]

    F_ij means the force between i and j, it's possitive, so we need a direction, 
    normalized(rj - ri) and i+ j - 
    """
    cR = np.array(cR)
    if not isinstance(cR, ndarray):
        raise TypeError(f'rs should be ndarray, not {type(cR).__name__}')
    r = np.linalg.norm(cR)

    beta = rc / sigma
    alpha = (27 / 4) * (beta ** 2 / (1 - beta ** 2) ** 3)
    return - epsilon * alpha * ((4 * rc**2 - 2 * sigma**2) * r**-3 \
                                             - (4 * rc**4 + 2 * sigma**2 * rc**2) * r**-5 
                                             + 4 * sigma**2 * rc**4 * r**-7) * cR / r

# update the particles' coordinates
def update_coord(cR, cV, cA, t, L) -> ndarray:
    t = float(t)
    cR = np.array(cR, dtype=float)
    cV = np.array(cV, dtype=float)
    cA = np.array(cA, dtype=float)
    nR = (cR + t * cV + t ** 2 * cA / 2) % L # a new method to preserve the pbc
    return nR

# update the particles' velocity
def update_velocity(cV, cA, nA, t) -> ndarray:
    return cV + t * (cA + nA) / 2

# generate a random maxewll distributed velocity
def maxwell_velocity(T: float, N: int, m: int) -> ndarray:
    """ using metropolis to generate a N dimensional velocity
        with maxwell distribution
    """
    rho = np.random.rand() * 2 * np.pi
    v0 = np.sqrt(2 * T) * np.array([np.cos(rho), np.sin(rho)])
    res = [v0]

    n = 0
    while n < m - 1:
        pv = pertubation(v0, T * 1e-2)
        e1 = np.linalg.norm(v0) ** 2
        e2 = np.linalg.norm(pv) ** 2
        if metropolis((e2 - e1) / 2, T):
            n += 1
            v0 = pv
            res.append(v0)
    return res

# calculate the potential energy of the ceter 
def center_potential(center, nearby_index, cell):
    u = 0
    for item in nearby_index: # all nearby cells
        try: 
            for particle in cell[str(item)]: # all particles in the cell
                u += potential(particle, center)
        except KeyError:
            continue 
    return u

def nearby(cell_index: str, size: int=999):
    """return nearby index of the cell"""
    key = cell_index.replace('\'', '\"')
    i, j = json.loads(key)
    nearby_index = [[i - 1, j], [i + 1, j], 
                    [i, j + 1], [i, j - 1], 
                    [i + 1, j + 1], [i + 1, j - 1], 
                    [i - 1, j + 1], [i - 1, j - 1]]
    pbc_index = []
    for item in nearby_index:
        i, j = item
        # periodical boundary condition
        if i < 0: i += size # this cannot be removed
        if i > size: i -= size
        if j < 0: j += size
        if j >= size: j -= size
        pbc_index.append([i, j]) 
    return pbc_index

def potential(particle: List[float], center: List[float]):
    """calcuate the interaction of two coordinates"""
    particle = np.array(particle)
    center = np.array(center)
    rc = 1.5 # cutoff r
    r = np.linalg.norm(particle - center)
    beta = rc / 1
    alpha = (27 / 4) * (beta ** 2 / (1 - beta ** 2) ** 3)
    if r <= rc:
        return - alpha * (1 / r **2 - 1) * (beta ** 2 / r ** 2 - 1) ** 2
    else:
        return 0


def pertubation(center: List[float], d: float=0.5):
    """return a pertubation of the center coordinates

    >>> center = [.13, .14]
    >>> pert_center = pertubation(center)
    >>> len(pert_center)
    2
    """
    # d = np.sqrt(L ** 2 (1 - rho) / N)
    if len(center) != 2:
        raise ValueError("center must be 2-dimensional")
    center = np.array(center)
    center[0] = center[0] + d * (np.random.rand() - .5)
    center[1] = center[1] + d * (np.random.rand() - .5)
    return center

def metropolis(du, T):
    """return a metropolis criterion"""
    if du < 0:
        return 1
    elif np.exp(- du / T) > np.random.rand():
        return 1
    else:
        return 0

def visualize(cell: dict):
    answer = []
    for key, value in cell.items():
        for item in value:
            answer.append(item)
    coordx = [item[0] for item in answer]
    coordy = [item[1] for item in answer]

    plt.scatter(coordx, coordy, s=1)
    plt.title("N = "+str(len(answer)))
    plt.show()

def timing(n, m):
    if n == int (m * 0.25): print("25%")
    if n == int (m * 0.50): print("50%")
    if n == int (m * 0.75): print("75%")
    if n == int (m * .99): print("99%")



# test how it works

if __name__ == '__main__':
    L = 17
    r = 0.5
    T = 1
    rho = 0.8442
    period = 0.1
    step = 100

    N = L * L * rho
    # d =  np.sqrt(L ** 2 (1 - rho) / N)
    coord, num_placed, cell= random_config(N, r, L)

    energy = monte_carol(cell, 0.7, 0, 1e5)
    # coord, velocity, energy = velocity_verlet(coord, maxwell_velocity(T, 2, N), cut_off_potential, period, step, L)

    # coord_old = np.array(coord, dtype=float)
    # plt.scatter(coord_old[:,0], coord_old[:,1], s=1)
    # plt.show()
    plt.plot(energy)
    plt.show()
    # plt.scatter(coord[:,0], coord[:,1], s=10)
    # plt.scatter(coord_old[:,0], coord_old[:,1], color='r', s=1)
    # print("N="+str(len(coord)))
    # plt.show()




    # print("total particles inserted : ", num_placed)
    # print("expected number of particles inserted : ", N)
    # # print(coord)
    # coordx = [item[0] for item in coord]
    # coordy = [item[1] for item in coord]

    # plt.scatter(coordx, coordy, s=1)
    # plt.title("N = "+str(num_placed))
    # plt.show()