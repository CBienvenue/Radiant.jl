# CP

The `CP` solver is the collision-probability method. It expands the volumetric and surface fluxes in Legendre polynomials of arbitrary order, assembles the collision, leakage and transmission probability matrices, and solves the resulting linear system directly (`"global"` mode) or by an interface-current sweep (`"sweeping"` mode). It is currently available in one-dimensional Cartesian geometry for the Boltzmann transport equation, with vacuum, reflective and surface-source boundaries.

## Structure
```@docs
Radiant.CP
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.CP,particle::Particle)
Radiant.set_solver_type(this::Radiant.CP,solver_type::String)
Radiant.set_legendre_order(this::Radiant.CP,volume_order::Int64,surface_order::Int64)
Radiant.set_surface_order(this::Radiant.CP,surface_order::Int64)
Radiant.set_mode(this::Radiant.CP,mode::String)
Radiant.set_acceleration(this::Radiant.CP,acceleration::String,parameter::Int64)
Radiant.set_convergence_criterion(this::Radiant.CP,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.CP,maximum_iteration::Int64)
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.CP)
Radiant.get_solver_type(this::Radiant.CP)
Radiant.get_legendre_order(this::Radiant.CP)
Radiant.get_surface_order(this::Radiant.CP)
Radiant.get_mode(this::Radiant.CP)
Radiant.get_schemes(this::Radiant.CP,geometry::Radiant.Geometry,isFC::Bool)
Radiant.get_is_full_coupling(this::Radiant.CP)
Radiant.get_acceleration(this::Radiant.CP)
Radiant.get_gmres_restart(this::Radiant.CP)
Radiant.get_anderson_depth(this::Radiant.CP)
Radiant.get_convergence_criterion(this::Radiant.CP)
Radiant.get_maximum_iteration(this::Radiant.CP)
```
