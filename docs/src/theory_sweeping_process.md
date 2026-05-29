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

## 1.6.4 Convergence Acceleration

The source iteration of Section 1.6.1 is a **fixed-point iteration**. Writing the within-group scattering source as the linear operator $\mathcal{S}$ and the sweep (inversion of $\mathcal{L}_{p,g}$) as $\mathcal{L}^{-1}$, one iteration reads

$$\Phi^{(k+1)} = \mathcal{L}^{-1}\mathcal{S}\,\Phi^{(k)} + \mathcal{L}^{-1}\!\left(Q^{\texttt{B,out-of-group}} + Q^{\texttt{ext}}\right) \;\equiv\; \mathcal{A}\,\Phi^{(k)} + b \,,$$

where the constant term $b$ gathers the out-of-group, external surface and (for CSD) incoming-energy sources. The iteration converges geometrically with the spectral radius $\rho(\mathcal{A})$, which approaches the within-group scattering ratio $c = \Sigma_{s,0}^{g\to g}/\Sigma_{t}^{g}$ in optically thick media. The number of iterations needed to reach a tolerance $\varepsilon$ then grows like $\ln\varepsilon / \ln c$, so as $c \to 1$ the unaccelerated iteration becomes prohibitively slow. Radiant offers four acceleration options through `set_acceleration`.

### Livolant (`"livolant"`)

Livolant's vector-extrapolation method [hebert2009applied](@cite) recombines the last three successive iterates $\Phi^{(k-2)}$, $\Phi^{(k-1)}$, $\Phi^{(k)}$. Defining the residuals

$$\Delta_{1} = \Phi^{(k-1)} - \Phi^{(k-2)} \,,\quad \Delta_{2} = \Phi^{(k)} - \Phi^{(k-1)} \,,$$

an accelerated iterate is constructed by minimizing $\|\Delta_{2} - \beta\Delta_{1}\|$ and forming the linear combination

$$\tilde{\Phi}^{(k)} = \Phi^{(k)} + \frac{\beta}{1-\beta}\,\Delta_{2} \,,$$

applied periodically every few source iterations. It is inexpensive and effective on moderately scattering-dominated problems.

### Anderson (`"anderson"`)

Anderson acceleration [anderson1965iterative, walker2011anderson](@cite) generalizes Livolant to a memory of depth $m$ (the optional `set_acceleration` parameter, default $3$). At iteration $k$ it forms the residual $f^{(k)} = \mathcal{A}\Phi^{(k)} + b - \Phi^{(k)}$ and computes mixing coefficients $\gamma$ by the least-squares problem

$$\gamma = \arg\min_{\gamma}\;\Big\| f^{(k)} - \textstyle\sum_{i=1}^{m}\gamma_{i}\,\big(f^{(k-m+i)} - f^{(k-m+i-1)}\big) \Big\|_{2} \,,$$

from which the next iterate is built as a damped combination of the stored iterates. Because it accelerates the exact fixed-point sequence, Anderson is a robust drop-in replacement that applies unchanged to every boundary condition and physics mode (BTE / BFP / CSD / Fokker–Planck).

### GMRES (`"gmres"`) and BiCGStab (`"bicgstab"`)

Instead of iterating the fixed point, the equivalent linear system

$$(\mathcal{I} - \mathcal{A})\,\Phi = b$$

is solved directly by a **matrix-free Krylov method** [saad1986gmres, vandervorst1992bicgstab, warsa2004krylov](@cite). Each matrix–vector product

$$(\mathcal{I} - \mathcal{A})\,v = v - \mathcal{L}^{-1}\mathcal{S}\,v$$

costs exactly one transport sweep, so the operator $\mathcal{A}$ is never assembled. Krylov methods converge in a number of sweeps that is essentially **independent of $\rho$**, which makes them dramatically more efficient than source iteration as $c \to 1$. The Krylov unknown includes the incoming boundary angular fluxes on reflective and periodic faces, so the acceleration is exact for every boundary condition. `GMRES` is restarted every `set_acceleration` parameter (default $30$) Krylov vectors and minimizes the residual over the whole subspace; `BiCGStab` uses short recurrences with a fixed, smaller memory footprint (two sweeps per iteration).

### Choosing a method

The benefit of acceleration is governed by the *within-group* scattering ratio $c$, not by the total scattering:

- **Electron BFP/CSD transport** moves most particles to a *lower* energy group (out-of-group source) or loses energy continuously through the CSD operator, both handled directly by the sweep. The within-group ratio $c$ is small, $\rho$ is small, and source iteration already converges in a handful of iterations — `"none"` or `"livolant"` is sufficient and the Krylov methods bring little gain.
- **Highly scattering, optically thick BTE problems** (e.g. neutral-particle transport in thick low-absorption media, and more generally 2D/3D problems where $\rho$ is larger) have $c \to 1$. Here Krylov/Anderson reduce the iteration count by one to two orders of magnitude.

As an illustration, for a one-group reflecting slab with isotropic scattering (so $\rho \approx c$), the number of in-group iterations to reach $\varepsilon = 10^{-7}$ is:

| $c\;(\approx\rho)$ | `none` | `livolant` | `anderson` | `gmres` | `bicgstab` |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.90  | 141   | 53    | 25 | 7 | 12 |
| 0.99  | 1 207 | 221   | 34 | 7 | 12 |
| 0.999 | 9 668 | 4 557 | 42 | 7 | 12 |

The GMRES sweep count stays essentially constant while source iteration scales as $1/(1-c)$. Near $\rho = 1$ the Krylov methods are also *more accurate*: the source-iteration stopping test on the flux change underestimates the true error by a factor $\sim 1/(1-\rho)$, whereas GMRES/BiCGStab drive the actual linear residual to the requested tolerance.

### Estimated spectral radius

Whatever the method, each in-group solve reports an **estimated spectral radius** $\hat{\rho}$, defined as the geometric-average relative-residual reduction per transport pass:

$$\hat{\rho} = \left(\frac{r_N}{r_0}\right)^{1/(N-1)} \,,$$

where $N$ is the number of passes performed and $r_0$, $r_N$ are the first and last relative residuals — $\lVert T(\Phi)-\Phi\rVert / \lVert\Phi\rVert$ for the fixed-point methods (`none`, `livolant`, `anderson`) and $\lVert b - (\mathcal{I}-\mathcal{A})\Phi\rVert / \lVert b\rVert$ for the Krylov methods (`gmres`, `bicgstab`). The estimate is returned as `NaN` when the group converges without iterating (for instance a group with no in-group source).

For plain source iteration ($\mathtt{none}$) this recovers the true spectral radius $\rho(\mathcal{A}) \approx c$ of the iteration operator. For the accelerated methods $\hat{\rho}$ instead measures the *effective* per-pass reduction, which is much smaller — illustrating the speed-up directly. For the one-group reflecting slab of the table above:

| $c\;(\approx\rho)$ | `none` | `livolant` | `anderson` | `gmres` | `bicgstab` |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.90  | 0.89 | 0.72 | 0.50 | 0.003 | 0.04 |
| 0.99  | 0.99 | 0.93 | 0.61 | 0.005 | 0.08 |
| 0.999 | 1.00 | 0.95 | 0.67 | 0.008 | 0.11 |

The estimate is printed at the end of each group solve and is available afterwards, per energy group, through the `get_spectral_radius` method of the `Computation_Unit` (Section 6).

## 1.6.5 Boundary Conditions in the Sweep

The boundary moments required at the start of a sweep along each direction are determined by the global boundary conditions of the geometry (Section 5):

- `"void"`: zero incoming flux.
- `"reflective"`: the incoming flux on a face is set equal to the outgoing flux of the mirror-symmetric direction, which is computed and stored during the same outer iteration.
- `"periodic"`: the incoming flux on a face is set equal to the outgoing flux on the opposite face of the same direction.

When a surface source is defined (Section 7.2), its angular distribution contributes additively to the incoming boundary moments.
