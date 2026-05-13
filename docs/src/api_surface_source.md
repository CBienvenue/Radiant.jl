# Surface_Source

## Structure
```@docs
Radiant.Surface_Source
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.Surface_Source,particle::Particle)
Radiant.set_intensity(this::Radiant.Surface_Source,intensity::Real)
Radiant.set_energy_group(this::Radiant.Surface_Source,energy_group::Int64)
Radiant.set_direction(this::Radiant.Surface_Source,direction::Vector{Float64})
Radiant.set_location(this::Radiant.Surface_Source,location::String)
Radiant.set_boundaries(this::Radiant.Surface_Source,axis::String,boundaries::Vector{Float64})
Radiant.set_legendre_order(this::Radiant.Surface_Source,legendre_order::Int64)
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.Surface_Source)
Radiant.get_normalization_factor(this::Radiant.Surface_Source)
Radiant.get_legendre_order(this::Radiant.Surface_Source)
```
