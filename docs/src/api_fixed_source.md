# Fixed_Sources

## Structure
```@docs
Radiant.Fixed_Sources
```

## Build
```@docs
Radiant.build(this::Radiant.Fixed_Sources)
```

## Setters
```@docs
Radiant.add_source(this::Radiant.Fixed_Sources,fixed_source::Union{Radiant.Surface_Source,Radiant.Volume_Source})
```

## Getters
```@docs
Radiant.get_source(this::Radiant.Fixed_Sources,particle::Particle)
Radiant.get_particles(this::Radiant.Fixed_Sources)
Radiant.get_normalization_factor(this::Radiant.Fixed_Sources)
```
