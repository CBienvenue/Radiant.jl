# Cross_Sections

## Structure
```@docs
Radiant.Cross_Sections
```

## Build and File I/O
```@docs
Radiant.build(this::Radiant.Cross_Sections)
Radiant.write(this::Radiant.Cross_Sections,file::String)
```

## Setters
```@docs
Radiant.set_source(this::Radiant.Cross_Sections,source::String)
Radiant.set_file(this::Radiant.Cross_Sections,file::String)
Radiant.set_materials(this::Radiant.Cross_Sections,materials::Union{Vector{Radiant.Material},Radiant.Material})
Radiant.set_particles(this::Radiant.Cross_Sections,particles::Union{Vector{Particle},Particle})
Radiant.set_energy(this::Radiant.Cross_Sections,energy::Vector{Float64})
Radiant.set_cutoff(this::Radiant.Cross_Sections,cutoff::Vector{Float64})
Radiant.set_number_of_groups(this::Radiant.Cross_Sections,number_of_groups::Vector{Int64})
Radiant.set_group_structure(this::Radiant.Cross_Sections,group_structure::Vector{T}) where T <: Real
Radiant.set_group_structure(this::Radiant.Cross_Sections,particle::Particle,group_structure::Vector{T}) where T <: Real
Radiant.set_group_structure(this::Radiant.Cross_Sections,type::String,Ng::Int64,E::Real,Ec::Real)
Radiant.set_group_structure(this::Radiant.Cross_Sections,particle::Particle,type::String,Ng::Int64,E::Real,Ec::Real)
Radiant.set_interactions(this::Radiant.Cross_Sections,interactions::Union{Vector{<:Radiant.Interaction},<:Radiant.Interaction})
Radiant.set_legendre_order(this::Radiant.Cross_Sections,legendre_order::Int64)
Radiant.set_absorption(this::Radiant.Cross_Sections,Σa::Vector{Float64})
Radiant.set_scattering(this::Radiant.Cross_Sections,Σs::Vector{Float64})
```

## Library Metadata Getters
```@docs
Radiant.get_file(this::Radiant.Cross_Sections)
Radiant.get_number_of_materials(this::Radiant.Cross_Sections)
Radiant.get_materials(this::Radiant.Cross_Sections)
Radiant.get_number_of_particles(this::Radiant.Cross_Sections)
Radiant.get_particles(this::Radiant.Cross_Sections)
Radiant.get_particle(this::Radiant.Cross_Sections,tag::String)
Radiant.get_material(this::Radiant.Cross_Sections,tag::String)
Radiant.get_densities(this::Radiant.Cross_Sections)
Radiant.get_energy(this::Radiant.Cross_Sections)
Radiant.get_cutoff(this::Radiant.Cross_Sections)
Radiant.get_legendre_order(this::Radiant.Cross_Sections)
Radiant.get_group_structure(this::Radiant.Cross_Sections)
Radiant.get_interactions(this::Radiant.Cross_Sections)
```

## Per-Particle Energy Structure
```@docs
Radiant.get_number_of_groups(this::Radiant.Cross_Sections)
Radiant.get_number_of_groups(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_energy_boundaries(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_energies(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_energy_width(this::Radiant.Cross_Sections,particle::Particle)
```

## Multigroup Data Getters
```@docs
Radiant.get_absorption(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_total(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_scattering(this::Radiant.Cross_Sections,particle_in::Particle,particle_out::Particle,legendre_order::Int64)
Radiant.get_stopping_powers(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_boundary_stopping_powers(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_momentum_transfer(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_energy_deposition(this::Radiant.Cross_Sections,particle::Particle)
Radiant.get_charge_deposition(this::Radiant.Cross_Sections,particle::Particle)
```
