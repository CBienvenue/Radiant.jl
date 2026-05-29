# SN

The `SN` struct describes the discrete-ordinates discretization method used for one particle. The alias `Discrete_Ordinates` is kept for backward compatibility.

## Structure
```@docs
Radiant.SN
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.SN,particle::Particle)
Radiant.set_solver_type(this::Radiant.SN,solver_type::String)
Radiant.set_quadrature(this::Radiant.SN,type::String,order::Int64,Qdims::Int64)
Radiant.set_legendre_order(this::Radiant.SN,legendre_order::Int64)
Radiant.set_angular_fokker_planck(this::Radiant.SN,angular_fokker_planck::String)
Radiant.set_angular_boltzmann(this::Radiant.SN,angular_boltzmann::String)
Radiant.set_convergence_criterion(this::Radiant.SN,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.SN,maximum_iteration::Int64)
Radiant.set_scheme(this::Radiant.SN,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_acceleration(this::Radiant.SN,acceleration::String,parameter::Int64)
Radiant.set_is_full_coupling(this::Radiant.SN,isFC::Bool)
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.SN)
Radiant.get_solver_type(this::Radiant.SN)
Radiant.get_quadrature_type(this::Radiant.SN)
Radiant.get_quadrature_order(this::Radiant.SN)
Radiant.get_quadrature_dimension(this::Radiant.SN)
Radiant.get_legendre_order(this::Radiant.SN)
Radiant.get_angular_boltzmann(this::Radiant.SN)
Radiant.get_angular_fokker_planck(this::Radiant.SN)
Radiant.get_convergence_criterion(this::Radiant.SN)
Radiant.get_maximum_iteration(this::Radiant.SN)
Radiant.get_acceleration(this::Radiant.SN)
Radiant.get_gmres_restart(this::Radiant.SN)
Radiant.get_anderson_depth(this::Radiant.SN)
Radiant.get_is_full_coupling(this::Radiant.SN)
Radiant.get_schemes(this::Radiant.SN,geometry::Radiant.Geometry,isFC::Bool)
```
