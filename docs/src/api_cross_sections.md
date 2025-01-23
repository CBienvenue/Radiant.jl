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
Radiant.set_energy(this::Radiant.Cross_Sections,energy::Real)
Radiant.set_cutoff(this::Radiant.Cross_Sections,cutoff::Real)
Radiant.set_number_of_groups(this::Radiant.Cross_Sections,number_of_groups::Vector{Int64})
Radiant.set_group_structure(this::Radiant.Cross_Sections,group_structure::Vector{String})
Radiant.set_interactions(this::Radiant.Cross_Sections,interactions::Vector{Radiant.Interaction})
Radiant.set_legendre_order(this::Radiant.Cross_Sections,legendre_order::Int64)
Radiant.build(this::Radiant.Cross_Sections)
Radiant.write(this::Radiant.Cross_Sections,file::String)
```