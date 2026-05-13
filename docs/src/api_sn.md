## Structure
```@docs
Radiant.SN
```

## Methods
```@docs
Radiant.set_particle(this::Radiant.SN,particle::Particle)
Radiant.set_solver_type(this::Radiant.SN,solver_type::String)
Radiant.set_quadrature(this::Radiant.SN,type::String,order::Int64)
Radiant.set_legendre_order(this::Radiant.SN,legendre_order::Int64)
Radiant.set_angular_fokker_planck(this::Radiant.SN,angular_fokker_planck::String)
Radiant.set_angular_boltzmann(this::Radiant.SN,angular_boltzmann::String)
Radiant.set_convergence_criterion(this::Radiant.SN,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.SN,maximum_iteration::Int64)
Radiant.set_scheme(this::Radiant.SN,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_acceleration(this::Radiant.SN,acceleration::String)
```