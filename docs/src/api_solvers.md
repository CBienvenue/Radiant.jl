# Solvers

A `Solvers` object groups one discretization method per particle and controls the coupling between them.

## Structure
```@docs
Radiant.Solvers
```

## Setters
```@docs
Radiant.add_solver(this::Radiant.Solvers,method::Union{Radiant.SN,Radiant.DPN,Radiant.GN})
Radiant.set_maximum_number_of_generations(this::Radiant.Solvers,maximum_number_of_generations::Int64)
Radiant.set_convergence_criterion(this::Radiant.Solvers,convergence_criterion::Real)
Radiant.set_convergence_type(this::Radiant.Solvers,convergence_type::String)
```

## Getters
```@docs
Radiant.get_method(this::Radiant.Solvers,particle::Particle)
Radiant.get_particles(this::Radiant.Solvers)
Radiant.get_number_of_particles(this::Radiant.Solvers)
Radiant.get_maximum_number_of_generations(this::Radiant.Solvers)
Radiant.get_convergence_criterion(this::Radiant.Solvers)
Radiant.get_convergence_type(this::Radiant.Solvers)
```
