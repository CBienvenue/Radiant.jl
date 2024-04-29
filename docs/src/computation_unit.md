## Structure
```@docs
Radiant.Computation_Unit
```

## Methods
```@docs
Radiant.set_cross_sections(this::Radiant.Computation_Unit,cross_sections::Radiant.Cross_Sections)
Radiant.set_geometry(this::Radiant.Computation_Unit,geometry::Radiant.Geometry)
Radiant.set_solvers(this::Radiant.Computation_Unit,methods::Radiant.Solvers)
Radiant.set_sources(this::Radiant.Computation_Unit,sources::Radiant.Fixed_Sources)
Radiant.run(this::Radiant.Computation_Unit)
Radiant.get_voxels_position(this::Radiant.Computation_Unit,axis::String)
Radiant.get_energies(this::Radiant.Computation_Unit,particle::String)
Radiant.get_flux(this::Radiant.Computation_Unit,particle::String)
Radiant.get_energy_deposition(this::Radiant.Computation_Unit,type::String)
Radiant.get_charge_deposition(this::Radiant.Computation_Unit,type::String)
```