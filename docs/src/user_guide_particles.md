# 3 Particles

In Radiant, particles are represented by the `Radiant.Particle` object. There are two primary ways to instantiate such an object:

**Using Reserved Particle Constructors:** Reserved constructors are designed for predefined particles like photons, electrons, or positrons. For example:
   
```julia
photon = Photon()
electron = Electron()
positron = Positron()
```

Each of these generates a `Radiant.Particle` object preloaded with the properties specific to that particle type.

**Using the General Constructor:** The general constructor allows the creation of custom particles. For instance:
   
```julia
photon = Particle()
```

This generates a Radiant.Particle object without any predefined properties.

These `Radiant.Particle` object has the following features:
- **Unique Identifiers:** Every instantiation of a `Radiant.Particle` object, whether using a reserved or general constructor, is assigned a unique ID. This ensures that each particle can be distinguished, even if their properties are identical.
- **Exclusive Reserved Particles:** Reserved particle types (e.g., photons, electrons) can only be instantiated through their specific constructors. This guarantees that these particles retain their predefined properties, enabling safe identification and treatment in computations.