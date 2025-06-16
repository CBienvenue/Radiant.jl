## Structure
```@docs
Radiant.Cross_Sections
```

## Methods
```@docs
Radiant.set_source(this::Radiant.Cross_Sections,source::String)
Radiant.set_file(this::Radiant.Cross_Sections,file::String)
Radiant.set_materials(this::Radiant.Cross_Sections,materials::Vector{Radiant.Material})
Radiant.set_particles(this::Radiant.Cross_Sections,particles::Vector{Particle})
Radiant.set_group_structure(this::Radiant.Cross_Sections,group_structure::Vector{T}) where T <: Real
Radiant.set_group_structure(this::Radiant.Cross_Sections,particle::Particle,group_structure::Vector{T}) where T <: Real
Radiant.set_group_structure(this::Radiant.Cross_Sections,type::String,Ng::Int64,E::Real,Ec::Real)
Radiant.set_group_structure(this::Radiant.Cross_Sections,particle::Particle,type::String,Ng::Int64,E::Real,Ec::Real)
Radiant.set_interactions(this::Radiant.Cross_Sections,interactions::Vector{Radiant.Interaction})
Radiant.set_legendre_order(this::Radiant.Cross_Sections,legendre_order::Int64)
Radiant.build(this::Radiant.Cross_Sections)
Radiant.write(this::Radiant.Cross_Sections,file::String)
```