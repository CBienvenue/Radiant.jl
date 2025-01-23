## Structure
```@docs
Radiant.Discrete_Ordinates
```

## Methods
```@docs
Radiant.set_particle(this::Radiant.Discrete_Ordinates,particle::Particle)
Radiant.set_solver_type(this::Radiant.Discrete_Ordinates,solver_type::String)
Radiant.set_quadrature(this::Radiant.Discrete_Ordinates,type::String,order::Int64)
Radiant.set_legendre_order(this::Radiant.Discrete_Ordinates,legendre_order::Int64)
Radiant.set_angular_fokker_planck(this::Radiant.Discrete_Ordinates,angular_fokker_planck::String)
Radiant.set_angular_boltzmann(this::Radiant.Discrete_Ordinates,angular_boltzmann::String)
Radiant.set_convergence_criterion(this::Radiant.Discrete_Ordinates,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.Discrete_Ordinates,maximum_iteration::Int64)
Radiant.set_scheme(this::Radiant.Discrete_Ordinates,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_acceleration(this::Radiant.Discrete_Ordinates,acceleration::String)
```