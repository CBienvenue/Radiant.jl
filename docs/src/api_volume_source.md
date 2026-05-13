# Volume_Source

## Structure
```@docs
Radiant.Volume_Source
```

## Setters
```@docs
Radiant.set_particle(this::Radiant.Volume_Source,particle::Particle)
Radiant.set_intensity(this::Radiant.Volume_Source,intensity::Real)
Radiant.set_energy_group(this::Radiant.Volume_Source,energy_group::Int64)
Radiant.set_boundaries(this::Radiant.Volume_Source,axis::String,boundaries::Vector{Float64})
```

## Getters
```@docs
Radiant.get_particle(this::Radiant.Volume_Source)
Radiant.get_normalization_factor(this::Radiant.Volume_Source)
```
