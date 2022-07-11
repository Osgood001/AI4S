# Molecular Simulation Note

[TOC]

## potentials

In general, the interaction potential can be expressed using LJ potential
$$
u_{LJ}(r_{ij})=4\epsilon \bigg [\bigg(\frac{\sigma}{r_{ij}}\bigg)^{12}- \bigg(\frac{\sigma}{r_{ij}}\bigg)^6\bigg ]
$$
But this is a long range potential which requires us to calculate the interaction of particles far away from each other, making it computationally costy.

Therefore, we use a more friendly potential with the following charaters:

1. the potential is repulsive when $r<\sigma$, and attractive within $\sigma < r <r_c$ range, where$r_c = k\sigma$, and $k=r_c/\sigma$ is not a very large number($k=2$ for atomic system and $k=1.2$ for colloidal system)
2.  $r_{\text{min}}\approx \sigma$
3. $\phi(r_{\text{min}})=\epsilon$

a simple choice would be
$$
\phi(r)=\epsilon \alpha \Big(\Big[ \frac{\sigma}{r}\Big]^2 -1\Big) \Big( \Big[ \frac{r_c}{r}\Big]^2 - 1\Big)^2
$$
if we set the unit that $\sigma = 1, \epsilon=1$, then the following can be simplified as
$$
\phi(r)=\alpha \Big(\Big[ \frac{1}{r}\Big]^2 -1\Big) \Big( \Big[ \frac{k}{r}\Big]^2 - 1\Big)^2
$$
we can obtain the expression of force
$$
F(r)=-\frac{d\phi(r)}{dr}=-\frac{2 a (k-r) (k+r) \left(k^2 \left(2 r^2-3\right)+r^2\right)}{r^7}
$$
given that $F(r_{\text{min}})=0$, we have
$$
r_{\text{min}}^2=3k^2/(2k^2+1)
$$
where $k=r_c/\sigma$

what's more, $\phi(r_{\text{min}})=\epsilon$, so
$$
\alpha = -\frac{27 k^2}{4 \left(k^2-1\right)^3}
$$
setting $k=2.5, 2, 1.2$, we have $\alpha = -0.29, -1, -114.11, r_{\text{min}}=1.18,1.15,1.06$ 

| Potential                | Lennard Jones $U^{LJ}(r)$ | ATOMIC $\phi^A(r)$ | COLLOIDAL $\phi^C(r)$ |
| ------------------------ | ------------------------- | ------------------ | --------------------- |
| Mininum $r_{\text{min}}$ | $1.12$                    | $1.15$             | $1.06$                |
| Cut-off $r_c$            | $2.5$                     | $2$                | $1.2$                 |
| Prefactor $\alpha$       | $N.A.$                    | $-1$               | $-114.11$             |

Compared from the table above, we see that LJ potential is similar to Atomic potential, while for colloidal, the prefactor is much larger, indicating a more narrow and deep potential well.(so colloidal are more likely to be stable than atomic configurations).

![image-20220704110536372](C:\Users\a1373\AppData\Roaming\Typora\typora-user-images\image-20220704110536372.png)

## Simulations

To simulate a particle system, one need to determine **the initial condition, the boundary condition and the simulation method**. 

For the **initial condition**, here we simulate a $L\times L$ box, with $N$ particles in it (Here we have inherited the unit by setting  $\sigma = 1, \epsilon=1$) and work with a canonical ensemble $NVT$.  We can start with a random configuration.

The **boundary condition** should be hard boundary in the $NVT$ case, but as we are dealing with $10^{23}$ particles in real system, but can only simulate $10^{3}$ particles, the results may be strongly related to the size (finite volumn effect).

The **simulation method** could be Monte Carol or Molecular Dynamics (the latter have many available algorithms), here we choose Monte Carol and Velocity-velet for our simulation.

### Initial condition

Randomly throw $N$ particles to $L\times L$ box can be a challenging task, so we apply the cell_list trick to accelate this process. The general idea would be

- SEPERATE: seperate the area into $N \times N$ boxes of length $d=1$ (diameter of the particle) and label them from $1$ to $N^2$ into an array

- SELECT: randomly select one box from the array

- THROW: try throwing a particle with its center in the box for several times(I throw $4000$ times)

  - if we find a valid position where it didn't collide with the nearby particles(which can be obtained by evaluating particles  only in the  nearby  $8$ boxes) then we remove this box from the array

  - if we didn't find a proper position, we will also remove this box from the array


- as the program proceeds, the array will have fewer element , so this algorithm converges really fast. it takes less than $1s$ to generate ($L=100, \rho = 0.7, N=7000$) scale without any error.

Note that it's not allowed to store $2$ particle center in one box (since there must be a distance of  $2$, but the largest length in the box would be $\sqrt{2}$ )



As for the lattice intial condtion, life gets easier by simply generate the initial coordinates. But to work with cell_list method (which will be useful in MC simulation), we need to convert the coord into a list recorded type.This is done by initialize the boxes and assign each particle to it.

 the general idea for lattice initiation would be:
```
- INPUT: N, L, r

- OUTPUT: coordinates, number_of_particles_placed, cell_list

  // generate lattice configuration

  - res = []

  - sqrtN = int (sqrt(N)) 

  - x = np.linspace(0, L, sqrtN)

  - y = product(x,x)

  - m = 0

  - for coord in y

    - res.append(coord)

    - m += 1

    - if (m == N):
      - break

  - dict = {}

  - for coord in res:

    - index = [int(res[0]), int(res[1])]

    - dict[str(index)] = res

  - return res, m, dict
```
### Boudary condition

During the simulation, it's possible that the coordinate of particles escape the box, a simple way is to do a modula on the coordinates with $L$, that is $\text{coord}_{\text{real}}=\text{coord}\mod{L}$ . In Python, $-8\mod{10}=2$, so this is exacty what we expect for PBC.

As for energy calculation, we need to refer to near by coordinates, this is done by setting  the index of box:

$[-1, n] \to [N - 1, n]$ (similar for other situations). 



### Simulation Method: Monte Carol

The general idea for MC method I implemented is shown below: 

```
- INPUT

  - res: the dict we obtained in the generation period, this data structure will be used generally as the configuration

  - T: the temperature of the monte carol simulation

  - u0: the energy of the configuration recorded by res (can be set to 0 since we care about relative energy)

  - m: step of valid MC needed, 100 as default

- n <- 0

- energy = [u0, ]

- while n < m:

  - reject <- 1

  - while  reject

    - choose a random cell from res

    - choose a random particle from cell

    - u1 <- calculate energy with surrounding particles

    - pertubate the particle (try to avoid collision)

    - u2 <- calculate energy with  pertubated surrounding particles

    - if metropolis(u2 - u1, T):

      - energy.append(u0 + u2 - u1)

      - reject <- 0

      - n += 1

      - update particle position 

- OUTPUT:

  - the energy sequence [u0, u1, u2, ... , um], where m is the number of valid moves in MC
```

### Simulation Method: Velocity-Verlet

The main part of this algorihtm is simple the iteration using 
$$
r_{k+1} =2 r_{k} - r_{k-1}+\tau^2 a_k\\
v_k = (r_{k+1} - r_{k - 1})/ (2\tau)
$$

```
- INPUT:

​    coord: coordinates of all particles

​    velocity: velocity of particles

​    potential: pairwise potential (refer to lj potential implemented before)

​    period: the length of one time step

​    step: the steps needed for the algorithm to stop

- OUTPUT:

​    final coordinates

​    final velocity

​    energy sequence

- for i in range(step):

  - cR = nR

  - cV = nV

  - cA = nA

  - nR = update_coord(cR, cV, cA, t, L)

  - nA = - force(nR) / mass

  - nV = update_velocity(cV, cA, nA, t)

  - energy.append(pair_potential_pbc(cR, cut_off_potential))

- return cR, cV, energy
```

But the PBC implementation is a bit different from the MC method:

- the energy and force calculation must be rewriten to calculate energies at boundary.

What's more, we need a Maxwell-Distribution velocity, which is implemented here using metropolis method.

```
using metropolis to generate a N dimensional velocity with maxwell distribution

INPUT: Temperature T, dimension N, particle number m

set initial velocity to v0
record v0
n < - 0
while n <m 
    pv <- perturbation(v0)
    if metropolis(energy difference, T)
        n +＝ 1
        v0 ＝ pv
        add velocity pv      
```

Another important detail is the thermo-stat, here we use Anderson thermo-stat.

## Results

**Task 2:**

Using Monte Carol, I obtained the following result(use random configuration) using $L=17, \rho = 0.8442, T=0.728$ and $L=40, \rho = 0.1, T=1.0$

![task1_first](D:\A0files\C1UCAS\D1study\E2Programming\Physics\Computational Physics\molecular liquids\results\task1_first.png)

![task1_second](D:\A0files\C1UCAS\D1study\E2Programming\Physics\Computational Physics\molecular liquids\results\task1_second.png)

We noticed that the average energy $U(T)=\frac{1}{T}\sum_t U(t)$  will converge as the simulation proceeds, a good criterium would be using $|U(T+\Delta T)- U(T)| < \delta$, where $\delta$ is the tolerance.



**Task 3:** **Finite size scaling**

Evaluate the $\left< U\right>$ and $\delta U^2$ for different system size $N$ .

Setting the size of the box $L=[10,20,30,40,50,60,70,80,90,100]$ and $T=0.728, \rho = 0.6$, we run MC simulation for 10000 times using random initial configuration.

<img src="C:\Users\a1373\AppData\Roaming\Typora\typora-user-images\image-20220704155118912.png" alt="image-20220704155118912" style="zoom: 33%;" /><img src="C:\Users\a1373\AppData\Roaming\Typora\typora-user-images\image-20220704154908649.png" alt="image-20220704154908649" style="zoom:33%;" />

The diagram is shown below:

![image-20220704155658442](C:\Users\a1373\AppData\Roaming\Typora\typora-user-images\image-20220704155658442.png)

which is almost linear under log scale, this indicate that $\delta U^2\propto L^m$, where $m\approx 4$ this indicate that
$$
\delta U \propto \sqrt{V}
$$
and the average energy is 

![image-20220704160438849](C:\Users\a1373\AppData\Roaming\Typora\typora-user-images\image-20220704160438849.png)

we can see that, as $L\to \infty$, $U/N$ approaches a constant value.

**Task 4: Radical distribution function**

In 2 dimensions, the rdf is $g(r)=\frac{dN(r)}{d(\pi r^2)}$, indicating the particle distribution of the system.

Here, the radical density function can be written as
$$
g(r)= \frac{2V}{N(N-1)Vr}\sum_{i=1}^N \sum_{j=1}^N \theta(r_{ij}-r)\theta(r + \Delta r - r_{ij})\\
\text{where~~~~} \theta(x)=
\begin{cases}
1, x> 0 \\
0,x\le 0
\end{cases}
$$
From book 4.51 we have
$$
U/N =  2\pi \rho \int_0^\infty dr r^2 u(r)g(r)
$$
where $\rho = \frac{\pi d^2 N}{4L^2}$

- $L=40, \rho = 0.1, T=1.0$

  the calculated average energy $u = U/N = - 0.17423$

  the measured average energy $u_m = U_m/N= - 1. 8389$

  we find that the relative error is $\delta u =2\frac{u_m - u}{u_m + u} = 165.38\%$
  
- $L=17, \rho = 0.8442, T=0.728$

  the calculated average energy $u = U/N = -3.092345$

  the measured average energy $u_m = U_m/N= - 2.19782$

  we find that the relative error is $\delta u =2\frac{u_m - u}{u_m + u} \approx 33.8\%$

- $L=10, \rho = 1.1, T=0.9$

  the calculated average energy $u = U/N = - 1.2389$

  the measured average energy $u_m = U_m/N= - 1.6267$

  we find that the relative error is $\delta u =2\frac{u_m - u}{u_m + u}\approx 27\%$

from the result, we see that when density increase, the relative error will be smaller, as the randomness and radical distribution wil statistically be more stable.

**Task 5: Equation of state**

From book 3.41,3.42 we have
$$
p =\frac{\rho}{\beta} + \frac{v_{ir}}{V}\\
v_{ir}=\frac{1}{2}\sum_i\sum_{j>i}f(r_{ij})r_{ij}
$$
here the coeffiecient is modified from $1/3$ to $1/2$ due to the system's two dimensionality.And the $f(r_{ij})$ is the force between $\text{i-th}$ and $\text{j-th}$ particle.

[Due to limited time to finish, there's no image here]

**Task 6: Free energy**

the phase diagram of the fliuid for the $\rho ,T$  with respect to $r_c$ 

[Due to limited time to finish, there's no image here]

for the LJ phase diagram is [taken from DOI:[10.1142/9789813209428_0007](http://dx.doi.org/10.1142/9789813209428_0007)]

![查看源图像](https://www.researchgate.net/profile/Athanasios_Fragkou/publication/5915231/figure/download/fig1/AS:281358258720775@1444092444739/Color-online-Phase-diagram-of-the-Lennard-Jones-system-from-14-The-simulated.png) 

 

