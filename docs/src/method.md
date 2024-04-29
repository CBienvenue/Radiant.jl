## Structure
```@docs
Radiant.Method
```

## Methods
```@docs
Radiant.set_particle(this::Radiant.Method,particle::String)
Radiant.set_solver_type(this::Radiant.Method,solver_type::String)
Radiant.set_quadrature(this::Radiant.Method,type::String,order::Int64)
Radiant.set_legendre_order(this::Radiant.Method,legendre_order::Int64)
Radiant.set_angular_fokker_planck(this::Radiant.Method,angular_fokker_planck::String)
Radiant.set_angular_boltzmann(this::Radiant.Method,angular_boltzmann::String)
Radiant.set_convergence_criterion(this::Radiant.Method,convergence_criterion::Float64)
Radiant.set_maximum_iteration(this::Radiant.Method,maximum_iteration::Int64)
Radiant.set_scheme(this::Radiant.Method,axis::String,scheme_type::String,scheme_order::Int64)
Radiant.set_acceleration(this::Radiant.Method,acceleration::String)
```