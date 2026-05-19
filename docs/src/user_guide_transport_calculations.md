# 8 Transport Calculations

The `Radiant.Computation_Unit` object brings together the cross-sections, geometry, solvers and sources defined in the previous chapters, executes the transport calculation, and exposes its results.

## 8.1 Assembling the Calculation

```julia
c = Computation_Unit()
c.set_cross_sections(cs)             # cs      :: Cross_Sections
c.set_geometry(geo)                  # geo     :: Geometry
c.set_solvers(solvers)               # solvers :: Solvers
c.set_sources(s)                     # s       :: Fixed_Sources
```

Each call assigns one of the four required components. The `Computation_Unit` ensures that they are consistent and complete before launching the calculation.

## 8.2 Running the Calculation

```julia
c.run()
```

When `run()` is called:

1. If the cross-section library has not yet been built (`cs.build()` not called), it is built automatically.
2. If the geometry has not yet been built, it is built automatically against the cross-sections.
3. The fixed sources are built automatically.
4. The transport calculation is executed and the resulting flux is stored inside the `Computation_Unit`.

Wrapping the call with Julia's `@time` macro provides timing information:

```julia
@time c.run()
```

## 8.3 Retrieving Results

After `run()` completes, several getter methods give access to the computed quantities.

### 8.3.1 Spatial and Energy Axes

```julia
x = c.get_voxels_position("x")       # Midpoint position of each voxel along x (cm)
y = c.get_voxels_position("y")       # In 2D and 3D
z = c.get_voxels_position("z")       # In 3D
E_electron = c.get_energies(electron) # Midpoint energy of each electron group (MeV)
E_photon   = c.get_energies(photon)
```

### 8.3.2 Energy and Charge Deposition

`get_energy_deposition` and `get_charge_deposition` return arrays whose shape matches the geometry voxelization. Without argument they return the total deposition over all particles; with a particle argument they return the contribution of that particle alone.

```julia
de_total    = c.get_energy_deposition()           # Total energy deposition
de_electron = c.get_energy_deposition(electron)   # Electron contribution only
dq_total    = c.get_charge_deposition()
dq_electron = c.get_charge_deposition(electron)
```

Energy deposition is given in MeV/g × cm² when the source intensity is unity, and charge deposition in cm²/g — they are normalized per source particle.

### 8.3.3 Angular-Flux Solution

The scalar flux per energy group and per voxel for a given particle is obtained with:

```julia
ϕ_electron = c.get_flux(electron)
ϕ_photon   = c.get_flux(photon)
```

The returned array has the structure `(Ng, Nx, Ny, Nz)` (with `Ny = Nz = 1` in lower dimensions).

## 8.4 Plotting Results

The results returned by `get_*` are standard Julia arrays and can be visualized with any plotting library (e.g. `Plots`, `PyPlot`, `Makie`). A typical 1D depth/dose plot is built as:

```julia
using PyPlot

x  = c.get_voxels_position("x")
de = c.get_energy_deposition()
dq = c.get_charge_deposition()

figure()
subplot(121); plot(x,de); xlabel("Depth (cm)"); ylabel("Energy deposition (MeV/g × cm²)")
subplot(122); plot(x,dq); xlabel("Depth (cm)"); ylabel("Charge deposition (cm²/g)")
```

## 8.5 Reusing a Calculation

The `Computation_Unit` does not store a deep copy of its inputs; modifying the cross-section, geometry, solver or source objects after a calculation will affect subsequent calls to `run()`. To run several variants of a calculation, instantiate a separate `Computation_Unit` per case or use deep copies of the input objects.

The cross-section library can be persisted to an FMAC-M file (see Section 4.4) and reused across calculations without rebuilding it:

```julia
cs.write("./library.txt")             # Save after a first run
# ... later ...
cs2 = Cross_Sections()
cs2.set_source("fmac-m")
cs2.set_file("./library.txt")
cs2.set_materials([water,al])
cs2.set_particles([electron,photon])
cs2.build()                           # Reload without recomputation
```

## 8.6 Summary of the Computation_Unit API

| Method                                       | Description                                                  |
|----------------------------------------------|--------------------------------------------------------------|
| `Computation_Unit()`                         | Constructor.                                                 |
| `set_cross_sections(cs)`                     | Assign the cross-sections library.                           |
| `set_geometry(geo)`                          | Assign the geometry.                                         |
| `set_solvers(solvers)`                       | Assign the solvers collection.                               |
| `set_sources(s)`                             | Assign the fixed sources collection.                         |
| `run()`                                      | Execute the transport calculation.                           |
| `get_voxels_position(axis)`                  | Voxel midpoint positions along an axis.                      |
| `get_energies(p)`                            | Group midpoint energies of particle `p`.                     |
| `get_energy_deposition([p])`                 | Total / per-particle energy deposition.                      |
| `get_charge_deposition([p])`                 | Total / per-particle charge deposition.                      |
| `get_flux(p)`                                | Scalar flux of particle `p` per group and per voxel.         |
