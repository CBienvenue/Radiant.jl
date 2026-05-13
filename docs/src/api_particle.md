# Particle

A `Particle` represents a particle type used in transport calculations. Reserved-particle constructors (`Photon`, `Electron`, ...) pre-load the mass and charge from physical constants; the generic `Particle(tag, mass, charge)` constructor creates a user-defined particle.

## Structure
```@docs
Radiant.Particle
```

## Setters
```@docs
Radiant.set_tag(this::Radiant.Particle,tag::String)
```

## Getters
```@docs
Radiant.get_tag(this::Radiant.Particle)
Radiant.get_mass(this::Radiant.Particle)
Radiant.get_charge(this::Radiant.Particle)
Radiant.get_type(this::Radiant.Particle)
```

## Type Predicates
```@docs
Radiant.is_photon(this::Radiant.Particle)
Radiant.is_electron(this::Radiant.Particle)
Radiant.is_positron(this::Radiant.Particle)
Radiant.is_proton(this::Radiant.Particle)
Radiant.is_antiproton(this::Radiant.Particle)
Radiant.is_alpha(this::Radiant.Particle)
Radiant.is_muon(this::Radiant.Particle)
Radiant.is_antimuon(this::Radiant.Particle)
```
