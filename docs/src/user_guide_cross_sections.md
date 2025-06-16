# 4 Cross-sections

A cross-sections library in Radiant is represented by the `Radiant.Cross_Sections` object. Once instantiated, its properties can be configured, and the necessary computations to format for particle transport can be performed.

## 4.1 Generating Cross-Sections from Physics Models

Radiant provides robust capabilities to produce accurate cross-sections for coupled or uncoupled transport of photons, electrons, and positrons. These cross-sections are derived from physics models and tabulated data. To generate them, several mandatory fields need to be specified:

```julia
cs = Cross_Sections()                       # Instantiate the cross-sections library
cs.set_source("physics-models")             # Specify cross-sections generation from physics models
cs.set_materials(material_list)             # Define the list of materials
cs.set_particles(particles_list)            # Define the list of particles
cs.set_group_structure("log",10,10.0,0.001) # Define the energy group structure (type, number of groups, midpoint energy of the highest energy group and cutoff energy). Specified that way, it defines the default discretization that applies for all particles.
cs.set_interactions(interactions_list)      # Define the list of physics interactions
cs.set_legendre_order(7)                    # Specify the Legendre moments order for scattering cross-sections
cs.build()                                  # Compute and build the transfer matrix for transport calculations
```

where
- `material_list` is a vector composed of `Radiant.Material` objects,
- `particles_list` is a vector composed of `Radiant.Particle` objects,
- `interactions_list` is a vector composed of `Radiant.Interaction` objects,
for example:

```julia
material_list = [water]
particles_list = [electron, photon]
interactions_list = [Inelastic_Leptons(),Elastic_Leptons(),Compton()]
```

**Custom Group Definitions:**
Radiant supports custom group confifgurations. It is sufficient then to directly gives the energy boundary vector, as follow:

```julia
cs.set_group_structure([0.1,0.5,1.0]) # or, equivalently cs.set_group_structure([1.0,0.5,0.1])
```

**Particle-Specific Group Definitions:**
Radiant supports particle-specific group configurations. The number of energy groups and the energy group structure can differ for each particle. For instance, if three particles are defined, their group settings might be configured as follows:

```julia
cs.set_group_structure(electron,"linear",10,10.0,0.1)
```

or 

```julia
cs.set_group_structure(electron,[0.1,0.5,1.0])
```

where `electron` is a `Radiant.Electron` object. This definition overwrite the default discretization if defined, which will be applied to all particles without explicit definition of group structure.

## 4.2 Generate Custom Cross-Sections

!!! note
    The capabilities to produce custom cross-sections are very limited at the moment.

Radiant provides some custom cross-sections capabilities. To generate them, several mandatory fields need to be specified:

```julia
cs = Cross_Sections()                     # Instantiate the cross-sections library
cs.set_source("custom")                   # Specify custom cross-sections
cs.set_materials([mat1,mat2,mat3,mat4])   # Define the list of materials
cs.set_particles([custom_particle])       # Define the list of particles
cs.set_absorption([50.0,5.0,0.0,0.1])     # Define the absorption cross-section per material
cs.set_scattering([0.0,0.0,0.0,0.9])      # Define the scattering (isotropic) cross-section per material
cs.build()                                # Compute and build the transfer matrix for transport calculations
```

These custom cross-sections capabilities are useful to test discretization schemes.

## 4.3 Read and Write Cross-Sections from FMAC-M files

Radiant is able to read and to write FMAC-M files. To read FMAC-M files and generate the formatted cross-sections for calculations:

```julia
cs = Cross_Sections()             # Instantiate the cross-sections library
cs.set_source("fmac-m")           # Specify the source of data as FMAC-M file
cs.set_file("./fmac_m_file.txt")  # Set the location and name of the FMAC-M file
cs.set_materials([water,al])      # Associate the data with material object.  
cs.build()
```

To write an FMAC-M file, it simply consists of

```julia
cs.write("./fmac_m_file.txt")     # Write the FMAC-M file. 
```