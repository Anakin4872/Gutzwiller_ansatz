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

$$\left.+ \sqrt{n+1}\sqrt{m} f_{n+1}^{(i)}\left(f_{m}^{*(i+1)} f_{m}^{(i+1)} + f_{m}^{*(i-1)} f_{m-1}^{(i-1)}\right)\right]$$

$$+ \frac{U}{2} f_n^{(i)}n(n-1)-\mu f_n^{(i)}n$$
    
 where the index $i$ runs over the $L$ lattice sites, $n$ and $m$ correspond to the internal degree of freedom (number of particles in this case), where $n,m={0, 1, 2}$ since the maximum allowed number of particles per site is two.

Later in this work, particles with spin-1 were considered, such that the index associated to the internal degrees of freedom is now associated to the possible magnetic projections of the particles at each site, 

$$\ket{\Psi}=\prod_i\sum_{m=-1}^{m=1} f_{m_i}^i \ket{m_i}$$


## **Use of the code**





