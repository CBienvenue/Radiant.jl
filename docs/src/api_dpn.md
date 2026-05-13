# DPN

The `DPN` solver is a moment-based ("double-Pn") angular discretization, used as an alternative to `SN`. It shares the same overall interface as `SN` but discretizes the angular variable on a polynomial basis instead of a quadrature.

## Structure
```@docs
Radiant.DPN
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.DPN,particle::Particle)
Radiant.set_solver_type(this::Radiant.DPN,solver_type::String)
Radiant.set_legendre_order(this::Radiant.DPN,legendre_order::Int64)
Radiant.set_polynomial_basis(this::Radiant.DPN,basis::String)
Radiant.set_angular_fokker_planck(this::Radiant.DPN,angular_fokker_planck::String)
Radiant.set_scheme(this::Radiant.DPN,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_is_full_coupling(this::Radiant.DPN,isFC::Bool)
Radiant.set_acceleration(this::Radiant.DPN,acceleration::String)
Radiant.set_convergence_criterion(this::Radiant.DPN,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.DPN,maximum_iteration::Int64)
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.DPN)
Radiant.get_solver_type(this::Radiant.DPN)
Radiant.get_legendre_order(this::Radiant.DPN)
Radiant.get_polynomial_basis(this::Radiant.DPN,Ndims::Int64)
Radiant.get_angular_fokker_planck(this::Radiant.DPN)
Radiant.get_schemes(this::Radiant.DPN,geometry::Radiant.Geometry,isFC::Bool)
Radiant.get_is_full_coupling(this::Radiant.DPN)
Radiant.get_acceleration(this::Radiant.DPN)
Radiant.get_convergence_criterion(this::Radiant.DPN)
Radiant.get_maximum_iteration(this::Radiant.DPN)
```
