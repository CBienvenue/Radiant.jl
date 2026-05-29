# GN

The `GN` solver is a moment-based angular discretization that subdivides the angular domain and uses a polynomial expansion on each subdivision. It generalizes the `DPN` solver and shares most of its interface.

## Structure
```@docs
Radiant.GN
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.GN,particle::Particle)
Radiant.set_solver_type(this::Radiant.GN,solver_type::String)
Radiant.set_legendre_order(this::Radiant.GN,legendre_order::Int64,legendre_order_local::Int64)
Radiant.set_subdivision(this::Radiant.GN,subdivision::Int64)
Radiant.set_polynomial_basis(this::Radiant.GN,basis::String)
Radiant.set_angular_fokker_planck(this::Radiant.GN,angular_fokker_planck::String)
Radiant.set_scheme(this::Radiant.GN,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_is_full_coupling(this::Radiant.GN,isFC::Bool)
Radiant.set_acceleration(this::Radiant.GN,acceleration::String,parameter::Int64)
Radiant.set_convergence_criterion(this::Radiant.GN,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.GN,maximum_iteration::Int64)
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.GN)
Radiant.get_solver_type(this::Radiant.GN)
Radiant.get_legendre_order(this::Radiant.GN)
Radiant.get_legendre_order_local(this::Radiant.GN)
Radiant.get_subdivision(this::Radiant.GN)
Radiant.get_polynomial_basis(this::Radiant.GN,Ndims::Int64)
Radiant.get_angular_fokker_planck(this::Radiant.GN)
Radiant.get_schemes(this::Radiant.GN,geometry::Radiant.Geometry,isFC::Bool)
Radiant.get_is_full_coupling(this::Radiant.GN)
Radiant.get_acceleration(this::Radiant.GN)
Radiant.get_gmres_restart(this::Radiant.GN)
Radiant.get_anderson_depth(this::Radiant.GN)
Radiant.get_convergence_criterion(this::Radiant.GN)
Radiant.get_maximum_iteration(this::Radiant.GN)
```
