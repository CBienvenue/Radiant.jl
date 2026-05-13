# Relaxation

The `Relaxation` interaction models both fluorescence and Auger cascades following inner-shell ionization. The convenience constructors `Fluorescence()` and `Auger()` activate only the radiative or non-radiative channels, respectively.

## Structure
```@docs
Radiant.Relaxation
```

## Setters
```@docs
Radiant.set_interaction_types(this::Radiant.Relaxation,interaction_types)
Radiant.set_minimum_probability(this::Radiant.Relaxation,ηmin::Real)
```
