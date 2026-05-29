# 6 Solvers

In Radiant, the discretization method used to transport a particle is described by a *solver* object. Radiant currently provides three solver types — `SN`, `DPN` and `GN` — each addressing a different angular-discretization paradigm. A `Solvers` object groups one solver per particle and controls the coupling between them.

## 6.1 The `Solvers` Object

A `Solvers` object is instantiated with no argument. Discretization methods are then added one by one with `add_solver`:

```julia
solvers = Solvers()
solvers.add_solver(sn_electron)
solvers.add_solver(sn_photon)
```

The order in which solvers are added is significant: it defines the order in which particles are processed during outer iterations of a coupled transport calculation.

### Coupled transport settings

Three additional methods control how the coupled outer iterations behave:

```julia
solvers.set_maximum_number_of_generations(20) # Maximum number of particle generations
solvers.set_convergence_criterion(1e-7)       # Outer-loop convergence criterion
solvers.set_convergence_type("flux")          # Quantity used to assess convergence
```

The convergence type can be:
- `"flux"` — converge the angular flux of each particle (default).
- `"energy-deposition"` — converge the total energy deposition.
- `"charge-deposition"` — converge the total charge deposition.

The maximum number of generations sets an upper bound on the number of secondary-particle cycles tracked when one particle produces another (e.g. electrons producing photons via bremsstrahlung, photons producing electrons via Compton, etc.).

## 6.2 Discrete Ordinates (`SN`) Solver

The `SN` solver is the workhorse of Radiant. It discretizes the angular variable on a quadrature set and the energy/spatial variables on user-specified schemes. The minimum information required to set up an `SN` solver is the particle, the solver type (i.e. which form of the transport equation), the angular quadrature, the Legendre order, and the spatial / energy schemes.

### 6.2.1 Form of the Transport Equation

```julia
sn.set_solver_type("BFP")
```

The available equations are:

| Value     | Equation                                                         |
|-----------|------------------------------------------------------------------|
| `"BTE"`   | Boltzmann transport equation                                     |
| `"BFP"`   | Boltzmann–Fokker–Planck equation                                 |
| `"BCSD"`  | Boltzmann with continuous slowing-down                           |
| `"FP"`    | Fokker–Planck only                                               |
| `"CSD"`   | Continuous slowing-down only                                     |
| `"BFP-EF"`| Boltzmann–Fokker–Planck without elastic                          |

### 6.2.2 Angular Quadrature

```julia
sn.set_quadrature("gauss-lobatto",8)
```

Available quadratures and their applicability:

| Quadrature                  | Applicability               |
|-----------------------------|-----------------------------|
| `"gauss-legendre"`          | 1D Cartesian only           |
| `"gauss-lobatto"`           | 1D Cartesian only           |
| `"carlson"`                 | 2D or 3D Cartesian          |
| `"gauss-legendre-chebychev"`| 2D or 3D Cartesian (product)|
| `"lebedev"`                 | 2D or 3D Cartesian          |

A third (optional) argument `Qdims` controls whether the quadrature spans the full sphere, the half sphere or only a line:

- `Qdims = 0` (default): infer from geometry,
- `Qdims = 1`: 1D quadrature over the direction cosine,
- `Qdims = 2`: half-sphere quadrature,
- `Qdims = 3`: full-sphere quadrature.

### 6.2.3 Legendre Expansion

The maximum order of the Legendre moments used to represent the scattering source is set with:

```julia
sn.set_legendre_order(7)
```

### 6.2.4 Angular Discretization of Boltzmann and Fokker–Planck Operators

```julia
sn.set_angular_boltzmann("galerkin")            # or "standard", "galerkin-m", "galerkin-d"
sn.set_angular_fokker_planck("finite-difference") # or "galerkin", "differential-quadrature"
```

The `"standard"` Boltzmann discretization is the classical discrete-ordinates approach. The Galerkin variants use a moment-to-discrete and discrete-to-moment couple of matrices: `"galerkin-d"` inverts the discrete-to-moment matrix (default), `"galerkin-m"` inverts the moment-to-discrete matrix, and `"galerkin"` is an alias for `"galerkin-d"`.

The Fokker–Planck operator can be discretized with:
- `"finite-difference"` — finite differences in angle (default),
- `"galerkin"` — moment-preserving Galerkin scheme,
- `"differential-quadrature"` — 1D Cartesian only.

### 6.2.5 Spatial and Energy Schemes

A scheme is set independently for each axis of the geometry and for the energy axis `"E"` (when continuous slowing-down is used):

```julia
sn.set_scheme("x","DG",2)   # Discontinuous Galerkin of order 2 along x
sn.set_scheme("E","DG",2)   # Discontinuous Galerkin of order 2 in energy
```

Available scheme types are:
- `"DD"` — diamond difference (any order),
- `"DG"` — discontinuous Galerkin (any order); also `"DG+"`, `"DG-"`,
- `"AWD"` — adaptive weighted difference (1st or 2nd order).

For multidimensional high-order schemes, the moments can be either fully coupled (e.g. `[00,10,01,11]`) or partially coupled (`[00,10,01]`). This is controlled with:

```julia
sn.set_is_full_coupling(true)   # default
```

### 6.2.6 In-Group Iteration Settings

```julia
sn.set_acceleration("livolant")    # in-group acceleration method
sn.set_convergence_criterion(1e-7) # default
sn.set_maximum_iteration(300)      # default
```

The available acceleration methods are:

| Value         | Method                                                        |
|---------------|---------------------------------------------------------------|
| `"none"`      | Plain source iteration (default).                             |
| `"livolant"`  | Livolant two-point vector extrapolation.                     |
| `"anderson"`  | Depth-$m$ Anderson acceleration of the fixed-point iteration. |
| `"gmres"`     | Matrix-free restarted GMRES on the equivalent linear system.  |
| `"bicgstab"`  | Matrix-free BiCGStab (lower memory than GMRES).               |

`"gmres"` and `"anderson"` accept an optional tuning parameter passed directly to
`set_acceleration` — the Krylov restart length for GMRES (default 30) and the memory depth for
Anderson (default 3):

```julia
sn.set_acceleration("gmres",50)    # GMRES restarted every 50 Krylov vectors
sn.set_acceleration("anderson",5)  # Anderson with memory depth 5
```

The convergence criterion and maximum iteration count are shared by all methods (for the Krylov solvers, the criterion is the relative-residual tolerance and the maximum iteration count bounds the number of transport sweeps).

!!! note
    The benefit of acceleration depends on the *within-group* scattering ratio. For electron BFP/CSD transport this ratio is small and source iteration already converges in a few iterations, so `"none"` or `"livolant"` is sufficient. The `"anderson"`, `"gmres"` and `"bicgstab"` methods pay off on highly scattering, optically thick BTE problems (and 2D/3D cases) where the spectral radius approaches one — see Section 1.6.4 of the Theory guide for a quantitative comparison.

### 6.2.7 Examples

1D BTE transport of photons with Galerkin quadrature:

```julia
sn = SN()
sn.set_particle(photon)
sn.set_solver_type("BTE")
sn.set_quadrature("gauss-lobatto",8)
sn.set_legendre_order(7)
sn.set_angular_boltzmann("galerkin")
sn.set_acceleration("livolant")
sn.set_convergence_criterion(4e-7)
sn.set_maximum_iteration(300)
sn.set_scheme("x","AWD",1)
```

3D BFP transport of electrons:

```julia
sn = SN()
sn.set_particle(electron)
sn.set_solver_type("BFP")
sn.set_quadrature("lebedev",8)
sn.set_legendre_order(15)
sn.set_angular_boltzmann("galerkin")
sn.set_angular_fokker_planck("finite-difference")
sn.set_acceleration("livolant")
sn.set_convergence_criterion(4e-7)
sn.set_maximum_iteration(300)
sn.set_scheme("E","DG",2)
sn.set_scheme("x","DD",1)
sn.set_scheme("y","DD",1)
sn.set_scheme("z","DD",1)
```

!!! note
    The `Discrete_Ordinates` name is a backward-compatible alias for `SN`. Legacy user code using `Discrete_Ordinates()` keeps working without modification.

## 6.3 `DPN` and `GN` Solvers

In addition to `SN`, Radiant provides two moment-based solver types: `DPN` (double-Pn) and `GN`. They share most of their interface with `SN` but replace the angular quadrature by a polynomial expansion on the angular variable.

The configuration steps are similar to `SN`:

```julia
dpn = DPN()
dpn.set_particle(electron)
dpn.set_solver_type("BFP")
dpn.set_legendre_order(7)
dpn.set_polynomial_basis("legendre")          # or "spherical-harmonics"
dpn.set_angular_fokker_planck("galerkin")
dpn.set_scheme("x","DG",2)
dpn.set_scheme("E","DG",2)
dpn.set_convergence_criterion(1e-7)
dpn.set_maximum_iteration(300)
dpn.set_acceleration("livolant")
dpn.set_is_full_coupling(true)
```

`GN` additionally exposes two methods specific to its angular treatment:

```julia
gn = GN()
gn.set_particle(electron)
gn.set_solver_type("BFP")
gn.set_polynomial_basis("legendre")           # 1D: Legendre-in-μ basis (default)
gn.set_legendre_order(16,2)                   # global order, local order
gn.set_subdivision(4)                         # number of angular sub-intervals
# ... other settings as for DPN
```

`set_legendre_order(L_global, L_local)` controls the global Legendre truncation order and a per-subdivision local order. `set_subdivision(n)` sets the number of angular subdivisions.

The polynomial basis available for `DPN` and `GN` are:
- `"legendre"` — a Legendre expansion in the polar cosine `μ` assuming azimuthal symmetry (default in 1D; for the `GN` solver this basis is available in 1D only),
- `"spherical-harmonics"` — a real spherical-harmonics expansion on the unit sphere (default in 2D and 3D; available in 1D, 2D and 3D).

For the `GN` solver in 1D, the `"legendre"` basis subdivides `μ ∈ [-1,1]` into `set_subdivision(n)` equal-width bands per hemisphere and carries a per-band local Legendre expansion of order `L_local`. It generalizes the `DPN` (double-Pn) discretization, which it reproduces with `set_subdivision(1)` and `L_local = L_global`.

!!! note
    The `DPN` solver supports only the `"galerkin"` discretization of the Fokker–Planck operator. The `GN` solver supports both `"galerkin"` and `"finite-difference"`; with the 1D `"legendre"` basis the latter uses a finite-volume discretization of the azimuthally-symmetric Laplace–Beltrami operator on the chain of `μ`-bands.

!!! note
    The `GN` solver supports the same in-group acceleration methods as `SN` — `"none"`, `"livolant"`, `"anderson"`, `"gmres"` and `"bicgstab"`, with the same optional `set_acceleration` tuning parameter (see Section 6.2.6). The `DPN` solver currently supports only `"none"` and `"livolant"`.

## 6.4 Mixing Solvers for Different Particles

A single calculation can mix solver types — for example using an `SN` solver for photons and a `GN` solver for electrons:

```julia
solvers = Solvers()
solvers.add_solver(sn_photon)
solvers.add_solver(gn_electron)
solvers.set_maximum_number_of_generations(10)
```

The `Solvers` object dispatches each particle to its associated method during coupled iterations.

## 6.5 Summary of the Solver API

### `Solvers`

| Method                                       | Description                                                |
|----------------------------------------------|------------------------------------------------------------|
| `Solvers()`                                  | Constructor.                                               |
| `add_solver(method)`                         | Register a particle-specific solver.                       |
| `set_maximum_number_of_generations(N)`       | Maximum number of generations in coupled transport.        |
| `set_convergence_criterion(ε)`               | Coupled-iteration convergence criterion.                   |
| `set_convergence_type(t)`                    | Convergence quantity (`"flux"`, etc.).                     |
| `get_method(p)`                              | Get the method associated with a particle.                 |
| `get_particles()`, `get_number_of_particles()`| Inspect the registered particles.                         |

### `SN`

| Method                                       | Description                                                |
|----------------------------------------------|------------------------------------------------------------|
| `SN()`                                       | Constructor.                                               |
| `set_particle(p)`                            | Particle the solver applies to.                            |
| `set_solver_type(t)`                         | `"BTE"`, `"BFP"`, `"BCSD"`, `"FP"`, `"CSD"`, `"BFP-EF"`.   |
| `set_quadrature(type,order[,Qdims])`         | Angular quadrature.                                        |
| `set_legendre_order(L)`                      | Legendre truncation order.                                 |
| `set_angular_boltzmann(t)`                   | Angular Boltzmann discretization.                          |
| `set_angular_fokker_planck(t)`               | Angular Fokker–Planck discretization.                      |
| `set_scheme(axis,type,order)`                | Spatial / energy scheme.                                   |
| `set_is_full_coupling(b)`                    | Toggle full coupling of high-order moments.                |
| `set_acceleration(t[,p])`                    | `"none"`, `"livolant"`, `"anderson"`, `"gmres"`, `"bicgstab"`; optional `p` = GMRES restart / Anderson depth. |
| `set_convergence_criterion(ε)`               | In-group convergence criterion.                            |
| `set_maximum_iteration(N)`                   | In-group maximum number of iterations.                     |

### `DPN` / `GN`

Same interface as `SN` for `set_particle`, `set_solver_type`, `set_legendre_order`, `set_scheme`, `set_acceleration`, `set_convergence_criterion`, `set_maximum_iteration`, `set_is_full_coupling`, plus:

| Method                                       | Description                                                |
|----------------------------------------------|------------------------------------------------------------|
| `set_polynomial_basis(b)`                    | `"legendre"` or `"spherical-harmonics"`.                   |
| `set_angular_fokker_planck("galerkin")`      | Currently the only available choice.                       |
| `set_subdivision(n)` *(GN only)*             | Number of angular subdivisions.                            |
| `set_legendre_order(L_g,L_l)` *(GN only)*    | Global and local Legendre orders.                          |
