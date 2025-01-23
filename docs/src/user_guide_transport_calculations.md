# 8 Transport Calculations

## 8.1 Executing Transport Calculations

A `Radiant.Computation_Unit` object has been defined in Radiant in order to launch transport calculation by putting together cross-sections, geometry, solvers and sources information. Once instantiated, the user can add these component and run the computation unit, for example:

```julia
c1 = Computation_Unit()
c1.set_cross_sections(cs)   # Set the cross-section object
c1.set_geometry(geo)        # Set the geometry object
c1.set_solvers(solvers)     # Set the solvers
c1.set_sources(s)           # Set the sources
c1.run()                    # Execute the transport calculations
``` 

This method ensure the completeness and consistency of the given data.

## 8.2 Analyzing Results

When transport calculations are done, the results are saved inside the `Radiant.Computation_Unit` object. It can be accessed using method, for example

```julia
x = c1.get_voxels_position("x")
energy_deposition = c1.get_energy_deposition()
charge_deposition = c1.get_charge_deposition()
```

to obtain the voxels mid-point along axis `x` and the corresponding energy or charge deposition. The flux solution for each particle can be accesses using

```julia
E_electron = c1.get_energies(electron)
flux_electron = c1.get_flux(electron)
```