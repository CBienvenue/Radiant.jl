# Computation_Unit

## Structure
```@docs
Radiant.Computation_Unit
```

## Setters
```@docs
Radiant.set_cross_sections(this::Radiant.Computation_Unit,cross_sections::Radiant.Cross_Sections)
Radiant.set_geometry(this::Radiant.Computation_Unit,geometry::Radiant.Geometry)
Radiant.set_solvers(this::Radiant.Computation_Unit,solvers::Radiant.Solvers)
Radiant.set_sources(this::Radiant.Computation_Unit,sources::Radiant.Fixed_Sources)
```

## Execution
```@docs
Radiant.run(this::Radiant.Computation_Unit)
```

## Result Getters
```@docs
Radiant.get_voxels_position(this::Radiant.Computation_Unit,axis::String)
Radiant.get_energies(this::Radiant.Computation_Unit,particle::Particle)
Radiant.get_flux(this::Radiant.Computation_Unit,particle::Particle)
Radiant.get_spectral_radius(this::Radiant.Computation_Unit,particle::Particle)
Radiant.get_energy_deposition(this::Radiant.Computation_Unit,particle::Particle)
Radiant.get_energy_deposition(this::Radiant.Computation_Unit)
Radiant.get_charge_deposition(this::Radiant.Computation_Unit,particle::Particle)
Radiant.get_charge_deposition(this::Radiant.Computation_Unit)
```
