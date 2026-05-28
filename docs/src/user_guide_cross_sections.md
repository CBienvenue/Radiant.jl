# 4 Cross-Sections

A cross-sections library in Radiant is represented by the `Radiant.Cross_Sections` object. Once instantiated, its properties are configured and the computations needed to format multigroup data for particle transport are performed with `cs.build()`.

Three sources of cross-sections are supported, and selected with `set_source(...)`:

| `set_source(...)` value | Description                                                              |
|-------------------------|--------------------------------------------------------------------------|
| `"physics-models"`      | Cross-sections are produced internally from physics models.              |
| `"fmac-m"`              | Cross-sections are read from a FMAC-M formatted file.                    |
| `"matxs"`               | Cross-sections are read from a MATXS (NJOY/CCCC) formatted file.          |
| `"custom"`              | User-supplied absorption and isotropic scattering cross-sections.        |

## 4.1 Generating Cross-Sections from Physics Models

Radiant produces accurate multigroup cross-sections for coupled or uncoupled transport of photons, electrons and positrons from a combination of physics models and tabulated atomic data.

A complete `physics-models` library requires the following inputs:

```julia
cs = Cross_Sections()
cs.set_source("physics-models")             # Cross-sections produced from physics models
cs.set_materials(material_list)             # List of Material objects
cs.set_particles(particle_list)             # List of Particle objects
cs.set_group_structure("log",20,10.0,0.001) # Default group structure for all particles
cs.set_interactions(interaction_list)       # List of Interaction objects
cs.set_legendre_order(7)                    # Truncation order for scattering Legendre moments
cs.build()                                  # Format the multigroup transfer matrices
```

where the inputs are:

- `material_list` — a vector of `Radiant.Material` objects (see Section 2),
- `particle_list` — a vector of `Radiant.Particle` objects (see Section 3),
- `interaction_list` — a vector of `Radiant.Interaction` objects (Section 4.4).

A typical configuration for coupled photon/electron transport in water reads:

```julia
material_list    = [water]
particle_list    = [electron, photon]
interaction_list = [Inelastic_Collision(),Elastic_Collision(),Bremsstrahlung(),
                    Compton(),Photoelectric(),Pair_Production(),
                    Annihilation(),Rayleigh(),Relaxation()]
```

### Energy group structure

Several forms of `set_group_structure` are available. The default discretization applies to every particle that does not have its own discretization defined.

**1. Pre-defined linear or logarithmic discretization:**
```julia
cs.set_group_structure("log",20,10.0,0.001)    # 20 logarithmically-spaced groups
cs.set_group_structure("linear",20,10.0,0.001) # 20 linearly-spaced groups
```
Arguments are `(type, number of groups, midpoint energy of the highest group [MeV], cutoff energy [MeV])`.

**2. Custom boundary vector:**
```julia
cs.set_group_structure([1.0,0.5,0.1])  # Two groups: [1.0,0.5] and [0.5,0.1]
cs.set_group_structure([0.1,0.5,1.0])  # Equivalent (any monotonic order is accepted)
```

**3. Particle-specific group structure:** Either of the previous forms can be applied to a specific particle by prepending the particle object. The particle-specific structure overrides the default for that particle only:
```julia
cs.set_group_structure(electron,"linear",80,10.0,0.001)
cs.set_group_structure(photon,[10.0,1.0,0.1,0.01,0.001])
```

The Legendre truncation order specifies the maximum order of the Legendre moments of the differential scattering cross-sections that will be computed and stored. The same Legendre order is used for all particles and all interactions.

## 4.2 The Interaction Library

Each physical process is represented by a subtype of `Interaction`. The available interactions are:

### Photons
- `Compton()` — Compton scattering. Models: `"klein-nishina"` (free electrons), `"waller-hartree"` (default, with incoherent scattering function).
- `Rayleigh()` — Coherent scattering on bound electrons.
- `Photoelectric()` — Photo-absorption. Models: `"epdl97"` (default, subshell dependent) or `"biggs_lighthill"`.
- `Pair_Production()` — e⁻/e⁺ pair creation. Angular distribution: `"sommerfield"` (default) or `"modified_dipole"`.

### Electrons and Positrons
- `Inelastic_Collision()` — Møller (e⁻) or Bhabha (e⁺) scattering. Density-effect corrections: `"fano"` (default), `"sternheimer"`, `"none"`. Subshell-dependent treatment and shell correction can be toggled. Solver target: `"BFP"` (default), `"FP"`, or `"BTE"`.
- `Elastic_Collision()` — Elastic scattering. Models: `"boschini"` (default, Mott) or `"rutherford"`. Optional Kawrakow and Seltzer corrections; extended transport correction.
- `Bremsstrahlung()` — Radiative energy loss. Photon angular distribution: `"sommerfield"` (default) or `"modified_dipole"`.
- `Annihilation()` — Positron–electron annihilation.

### Atomic Relaxation
- `Relaxation()` — Both fluorescence and Auger cascades following inner-shell ionization.
- `Fluorescence()` — Convenience constructor; only the fluorescence channels of `Relaxation()`.
- `Auger()` — Convenience constructor; only the Auger-electron channels of `Relaxation()`.

### Customizing Interactions

Each interaction object exposes setter methods to tailor the model used. Two examples:

```julia
elastic = Elastic_Collision()
elastic.set_model("boschini",true,true)  # Boschini Mott with Kawrakow and Seltzer corrections
elastic.set_transport_correction(true)   # Apply extended transport correction
elastic.set_solver("BFP")                # Split between Boltzmann and Fokker-Planck operators

inelastic = Inelastic_Collision()
inelastic.set_density_correction("sternheimer")
inelastic.set_is_subshells_dependant(true)
inelastic.set_is_shell_correction(true)
inelastic.set_scattering_model("BFP")
```

The active production / scattering / absorption channels of an interaction are controlled with `set_interaction_types(dict)`, where `dict` maps a `(incoming_particle_type, outgoing_particle_type)` tuple to a list of single-letter tags:

- `"S"` — scattered (degraded) particle of the same type as the incident one.
- `"P"` — produced (secondary) particle.
- `"A"` — absorbed incident particle.

For example, to disable knock-on electron production in inelastic collisions while keeping scattering:

```julia
inelastic.set_interaction_types( Dict((Electron,Electron) => ["S"]) )
```

A more complete description of each interaction can be found in the Theory section and in the API reference.

## 4.3 Generating Custom Cross-Sections

!!! note
    Custom cross-sections are intended primarily to test discretization schemes; their capabilities are deliberately limited.

The `"custom"` source bypasses the physics models and lets the user specify absorption and isotropic scattering cross-sections per material:

```julia
cs = Cross_Sections()
cs.set_source("custom")
cs.set_materials([mat1,mat2,mat3,mat4])
cs.set_particles([custom_particle])
cs.set_absorption([50.0,5.0,0.0,0.1])    # Σa per material [cm⁻¹]
cs.set_scattering([0.0,0.0,0.0,0.9])     # Σs per material [cm⁻¹]
cs.build()
```

## 4.4 Reading and Writing FMAC-M Files

Radiant can both read and write FMAC-M formatted files. The FMAC-M file stores a fully-built multigroup library and can therefore be reused between calculations without recomputing the cross-sections.

To read an FMAC-M file:

```julia
cs = Cross_Sections()
cs.set_source("fmac-m")
cs.set_file("./fmac_m_file.txt")
cs.set_materials([water,al])              # Must match the materials stored in the file
cs.set_particles([electron,photon])       # Must match the particles stored in the file
cs.build()
```

To write the current cross-section library to disk:

```julia
cs.write("./fmac_m_file.txt")
```

`cs.write(...)` will call `cs.build()` automatically if the library has not yet been built. The
output format is selected with the optional second argument (`"fmac-m"` by default):

```julia
cs.write("./fmac_m_file.txt")            # FMAC-M (default)
cs.write("./library.matxs","matxs")      # MATXS, BCD/ASCII encoding
cs.write("./library.matxs","matxs",true) # MATXS, binary encoding
```

## 4.5 Reading and Writing MATXS Files

Radiant can also read and write the [MATXS](@ref "The MATXS File Format") format — the
generalized multigroup interface file based on the CCCC conventions and produced by NJOY's
`matxsr` module. As with FMAC-M, a MATXS file stores a fully-built multigroup library that can
be reused between calculations.

To read a MATXS file (the BCD/ASCII and binary encodings are auto-detected):

```julia
cs = Cross_Sections()
cs.set_source("matxs")
cs.set_file("./library.matxs")
cs.set_materials([water,al])              # Must match the materials stored in the file (in order)
cs.set_particles([electron,photon])       # Must match the particles stored in the file
cs.build()
```

To write the current library as a MATXS file:

```julia
cs.write("./library.matxs","matxs")       # BCD/ASCII encoding
cs.write("./library.matxs","matxs",true)  # binary encoding
```

The record structure, the reaction-name keywords (`*tot0`, `*heat`, `*scat`, `*edep`, …) and
the particle codes used by Radiant are detailed in [The MATXS File Format](@ref).

## 4.6 Retrieving Cross-Section Data

After `build()` has been called, the multigroup data can be inspected through the getter methods of `Cross_Sections`. They return arrays indexed by group (and material, where applicable):

```julia
Σt  = cs.get_total(electron)                # Total cross-sections [g×Ng, Nmat]
Σa  = cs.get_absorption(electron)           # Absorption cross-sections
Σs  = cs.get_scattering(electron,photon,7)  # Scattering matrix up to Legendre order 7
S   = cs.get_stopping_powers(electron)      # Restricted stopping powers
Sb  = cs.get_boundary_stopping_powers(electron)
T   = cs.get_momentum_transfer(electron)    # Momentum-transfer cross-sections
Σe  = cs.get_energy_deposition(electron)    # Energy-deposition cross-sections
Σc  = cs.get_charge_deposition(electron)    # Charge-deposition cross-sections
```

The energy group midpoints, boundaries and widths for a given particle can be obtained with:

```julia
E   = cs.get_energies(electron)             # Midpoint energy of each group
Eb  = cs.get_energy_boundaries(electron)    # Group boundaries
ΔE  = cs.get_energy_width(electron)         # Group widths
```

Particles and materials stored in a `Cross_Sections` library (especially relevant when loading from disk) can be retrieved by tag or, for particles, by abstract type:

```julia
e = cs.get_particle("electron")  # By tag
e = cs.get_particle(Electron)    # By type
w = cs.get_material("water")
```

## 4.7 Summary of the Cross_Sections API

| Method                                            | Description                                                  |
|---------------------------------------------------|--------------------------------------------------------------|
| `Cross_Sections()`                                | Constructor.                                                 |
| `set_source(s)`                                   | Source: `"physics-models"`, `"fmac-m"`, `"matxs"`, `"custom"`. |
| `set_file(path)`                                  | File path (used by `"fmac-m"` and `"matxs"`).                |
| `set_materials(list)`                             | Set the list of materials.                                   |
| `set_particles(list)`                             | Set the list of particles.                                   |
| `set_interactions(list)`                          | Set the list of interaction models.                          |
| `set_legendre_order(L)`                           | Maximum Legendre order for differential cross-sections.      |
| `set_group_structure(...)`                        | Set the energy group structure (multiple forms).             |
| `set_energy(E)`, `set_cutoff(Ec)`                 | Set per-particle highest energy / cutoff energy.             |
| `set_number_of_groups(Ng)`                        | Set per-particle number of energy groups.                    |
| `set_absorption(Σa)`, `set_scattering(Σs)`        | Set custom cross-sections.                                   |
| `build()`                                         | Build the multigroup transfer matrices.                      |
| `write(path[,format[,binary]])`                   | Write the library; `format` is `"fmac-m"` (default) or `"matxs"`. |
| `get_energies(p)`, `get_energy_boundaries(p)`, `get_energy_width(p)` | Query the energy structure for particle `p`.              |
| `get_total(p)`, `get_absorption(p)`               | Get total / absorption cross-sections.                       |
| `get_scattering(pin,pout,L)`                      | Legendre moments of the scattering matrix.                   |
| `get_stopping_powers(p)`, `get_boundary_stopping_powers(p)` | Stopping powers.                                  |
| `get_momentum_transfer(p)`                        | Momentum-transfer cross-sections.                            |
| `get_energy_deposition(p)`, `get_charge_deposition(p)` | Energy / charge deposition cross-sections.              |
| `get_particle(tag_or_type)`, `get_material(tag)`  | Retrieve a stored particle or material.                      |
