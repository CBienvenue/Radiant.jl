# 1.6 Sweeping Process

Once the BFP equation has been discretized in energy (Section 1.2), in angle (Section 1.3) and in space (Section 1.4), and once a closure relation has been chosen (Section 1.5), the resulting linear system can be solved on a voxel-by-voxel basis along each discrete direction $\mathbf{\Omega}_{n}$. This direction-by-direction, voxel-by-voxel resolution is the **sweep**.

## 1.6.1 The Source-Iteration Algorithm

The multigroup BFP equation can be written symbolically, for a particle $p$ and an energy group $g$, as

$$\mathcal{L}_{p,g}\,\Phi_{p,g} = Q^{\texttt{B,within-group}}_{p,g}[\Phi_{p,g}] + Q^{\texttt{B,out-of-group}}_{p,g}[\Phi_{p,g'\neq g}] + Q^{\texttt{ext}}_{p,g} \,,$$

where $\mathcal{L}_{p,g}$ contains the streaming, total cross-section and angular/continuous-slowing-down operators. Because the within-group scattering source depends on the unknown $\Phi_{p,g}$ itself, the equation is solved iteratively (source iteration):

$$\mathcal{L}_{p,g}\,\Phi^{(k+1)}_{p,g} = Q^{\texttt{B,within-group}}_{p,g}\!\left[\Phi^{(k)}_{p,g}\right] + Q^{\texttt{B,out-of-group}}_{p,g}\!\left[\Phi_{p,g'\neq g}\right] + Q^{\texttt{ext}}_{p,g} \,.$$

At each iteration, the right-hand side is known and the operator $\mathcal{L}_{p,g}$ is inverted by the sweep described below. Iterations continue until the in-group convergence criterion (`set_convergence_criterion`) is met, up to `set_maximum_iteration`. When Livolant acceleration is enabled (`set_acceleration("livolant")`), the iterates are recombined as described in Section 1.6.4.

## 1.6.2 Energy Group and Particle Ordering

Coupled transport of $|P|$ particles with $G$ energy groups is solved by nesting three loops:

- **Outer loop**: particle generations. Each generation transports every particle from the highest-energy group to the lowest, taking into account the secondary particles produced in previous generations. Iteration stops when no new generation produces fluxes above the convergence threshold (Section 6.1).
- **Particle loop**: in the order given to `Solvers.add_solver(...)`, each particle is transported group by group. The scattering source from other particles is taken from the previous generation.
- **Group loop**: groups are processed from $g=1$ (highest energy) to $g=G$ (lowest energy). For each group, the source-iteration algorithm of Section 1.6.1 is applied.

The CSD operator (when active) introduces a unidirectional coupling between consecutive groups: $\Phi_{p,g-1/2}$ feeds $\Phi_{p,g}$. This is consistent with the highest-to-lowest group ordering.

## 1.6.3 Direction-by-Direction Sweep

For a given direction $\mathbf{\Omega}_{n} = (\mu_{n},\eta_{n},\xi_{n})$, the streaming operator dictates the order in which voxels are visited:

- $\mu_{n} > 0$: sweep $x$ from $i=1$ to $i=N_{x}$; otherwise from $i=N_{x}$ to $i=1$.
- $\eta_{n} > 0$: sweep $y$ from $j=1$ to $j=N_{y}$; otherwise reversed.
- $\xi_{n} > 0$: sweep $z$ from $k=1$ to $k=N_{z}$; otherwise reversed.

For each voxel reached in this order, the closure relation expresses the outgoing boundary moments in terms of the incoming boundary moments (already known from upstream voxels or boundary conditions) and the volume moments. The discrete BFP system inside the voxel is then a small dense linear system in the unknowns

$$\left\{ \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})},\; \Phi_{p,g,n,i^{\pm},j,k}^{(\cdot)},\; \Phi_{p,g,n,i,j^{\pm},k}^{(\cdot)},\; \Phi_{p,g,n,i,j,k^{\pm}}^{(\cdot)},\; \Phi_{p,g^{\pm},n,i,j,k}^{(\cdot)} \right\} \,,$$

of size at most $\mathcal{O}_{x}\mathcal{O}_{y}\mathcal{O}_{z}\mathcal{O}_{\texttt{E}}$ (when high-order moments are fully coupled, see Section 1.5.4). This system is solved by direct inversion; the outgoing moments are then stored to feed the next voxel. The procedure is carried out independently for each direction $\mathbf{\Omega}_{n}$, which is the basis of the strong parallelism of discrete-ordinates methods.

## 1.6.4 Livolant Acceleration

Source iteration converges slowly when the scattering ratio approaches unity. Radiant implements Livolant's vector-extrapolation method [hebert2009applied](@cite), which accelerates convergence by recombining the last three successive flux iterates $\Phi^{(k-2)}$, $\Phi^{(k-1)}$, $\Phi^{(k)}$. Defining the residuals

$$\Delta_{1} = \Phi^{(k-1)} - \Phi^{(k-2)} \,,\quad \Delta_{2} = \Phi^{(k)} - \Phi^{(k-1)} \,,$$

an accelerated iterate is constructed by minimizing $\|\Delta_{2} - \beta\Delta_{1}\|$ in a chosen norm and forming the linear combination

$$\tilde{\Phi}^{(k)} = \Phi^{(k)} + \frac{\beta}{1-\beta}\,\Delta_{2} \,.$$

The acceleration is applied periodically every few source iterations once the residuals have stabilized. It is particularly effective on scattering-dominated problems (e.g. photons in low-Z materials) and on CSD-dominated electron transport.

## 1.6.5 Boundary Conditions in the Sweep

The boundary moments required at the start of a sweep along each direction are determined by the global boundary conditions of the geometry (Section 5):

- `"void"`: zero incoming flux.
- `"reflective"`: the incoming flux on a face is set equal to the outgoing flux of the mirror-symmetric direction, which is computed and stored during the same outer iteration.
- `"periodic"`: the incoming flux on a face is set equal to the outgoing flux on the opposite face of the same direction.

When a surface source is defined (Section 7.2), its angular distribution contributes additively to the incoming boundary moments.
