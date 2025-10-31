# **Gutzwiller Ansatz for the Bose-Hubbard model**

The Gutzwiller wave function is constructed by applying a Gutzwiller projection operator to a
product state of single-particle orbitals, which effectively suppresses double occupancy of elec
tronic states. The projection operator acts on each lattice site and reduces the probability of twoelectrons occupying the same site due to strong electron-electron repulsion.

 In its most basic form, this method relies on approximating the many-body wave function by multiplying individual site contributions:
 
 $$\ket{\Psi}=\prod_i\sum_{n=0}^{n_{max}} f_{n_i}^i \ket{n_i}$$
 
where $\ket{n_i}$ denotes the Fock state of $n$ atoms in the i-th lattice site, $n_{max}$ is a system size-independent cut off in the number of atoms per site, and $f_{n}^{i}$ corresponds to the amplitude of having $n$ atoms in the i-th lattice site, which in this case corresponds to complex time-dependent coefficients. The amplitudes are normalized to $\sum_n|f_{n_i}^i|^2=1$. In the dynamical case, these amplitudes can be found by considering the variational principle, minimizing of the action of the system giving by $S=\int dt \mathcal{L}$. The lagrangian associated to the Hamiltonian can be written as:

$$\mathcal{L}= \frac{i\hbar}{2}\left(\braket{\Psi|\dot{\Psi}}-\braket{\dot{\Psi}|\Psi} \right) -\bra{\Psi}\hat{H}\ket{\Psi}$$

The dynamical equations can be retrieved from the Euler-Lagrange formulation,

$$\frac{d}{dt}\left(\frac{d\mathcal{L}}{d\dot{f}_{n_i}^{*i}}\right)=\frac{d\mathcal{L}}{df_{n_i}^{*i}}$$

After implementing the aforementioned steps and performing the corresponding algebraic procedures for the spinless Bose-Hubbard model with nearest neighbor interaction, the dynamics can be characterized by the following set of coupled differential equations,

$$i\hbar\frac{d}{dt} f_n^{(i)}=-t\sum_m \left[\sqrt{n} \sqrt{m} f_{n-1}^{(i)} \left(f_{m-1}^{*(i+1)} f_{m}^{(i+1)}+ f_{m-1}^{*(i-1)} f_{m}^{(i-1)}\right)\right.$$

$$+ \sqrt{n+1}\sqrt{m} f_{n+1}^{(i)}\left(f_{m}^{*(i+1)} f_{m}^{(i+1)} + f_{m}^{*(i-1)} f_{m-1}^{(i-1)}\right)$$

$$+ \frac{U}{2} f_n^{(i)}n(n-1)-\mu f_n^{(i)}n$$
    
 where the index $i$ runs over the $L$ lattice sites, $n$ and $m$ correspond to the internal degree of freedom (number of particles in this case), where $n,m={0, 1, 2}$ since the maximum allowed number of particles per site is two.

Later in this work, particles with spin-1 were considered, such that the index associated to the internal degrees of freedom is now associated to the possible magnetic projections of the particles at each site, 

$$\ket{\Psi}=\prod_i\sum_{m=-1}^{m=1} f_{m_i}^i \ket{m_i}$$


# **Use of the code**
## **For BHspinless2.C**
This program numerically solves the spinless Bose–Hubbard model using the Gutzwiller mean-field ansatz. The implementation is written in C++ and relies on the Armadillo linear algebra library and OpenMP for parallel computation. The goal is to compute observables such as the average particle number and number fluctuations across a range of hopping amplitudes and chemical potentials.

The code evolves the local Gutzwiller wavefunction coefficients $f_i(n)$ using imaginary-time evolution and the Runge–Kutta 4th-order (RK4) method. Each site is represented by a complex matrix of coefficients, and self-consistency is achieved through iterative normalization of the wavefunction at every site.

### Key model parameters:
- $\mu$: Chemical potential (scanned from MUmin to MUmax)
- t: Hopping amplitude (scanned from Tmin to Tmax)
- V: Density–density interaction between neighboring sites
- n: Maximum occupation number per site
- s: Number of lattice sites
- h: RK4 time step
- jmax: Maximum number of imaginary-time iterations

The output consists of the average number of particles ⟨n⟩ and its standard deviation $Δn^2$ for each pair of parameters ($\mu$, t).

**Notes**

- The code uses periodic boundary conditions and includes nearest-neighbor density interactions through the function `ITEpV2`.
- Parallelization (#pragma omp parallel for) speeds up the computation of different hopping amplitudes.
- All parameters can be easily modified at the top of the `main()` function.
- The physical constants are dimensionless for simplicity.


## **For BilinearBiquadraticHeisenberg.cpp**
This program numerically solves the spin-1 Bilinear–Biquadratic Heisenberg model within the Gutzwiller mean-field framework. It is implemented in C++ using the Boost uBLAS linear algebra library and OpenMP for potential parallelization. The solver computes the ground-state energy and local spin-state populations across a range of Hamiltonian parameters.

The model Hamiltonian considered is:

$$ H = \sum_{\langle i,j \rangle} \left[\cos(\theta)\,(\mathbf{S}_i \cdot \mathbf{S}_j) +\sin(\theta)\,(\mathbf{S}_i \cdot \mathbf{S}_j)^2\right]+D \sum_i (S_i^z)^2$$

The wavefunction at each site is represented as:

$$\ket{\psi_i} = f_{i,-1}\ket{m=-1} + f_{i,0}\ket{m=0} + f_{i,+1}\ket{m=+1}, \qquad \sum_{m=-1}^{+1} \lvert f_{i,m}\rvert^2 = 1$$

The program scans over ranges of the Hamiltonian parameters:

- θ: interaction angle, in radians (set via th_)
- D: single-ion anisotropy strength

For each pair of values (θ, D):

1. The code initializes an equal-weight spin state $f_i(m)=1/ \sqrt{3}$
2. The system evolves under imaginary time using `RK4()` and the differential operator `ITEp()`.
3. The energy is evaluated iteratively with `Energy()` until convergence (difference ≤ 10⁻⁵).
4. Results are stored in `.dat` files containing the energy trajectory and final state populations.
