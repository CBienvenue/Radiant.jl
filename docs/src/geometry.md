## Structure
```@docs
Radiant.Geometry
```

## Methods
```@docs
Radiant.set_type(this::Radiant.Geometry,type::String)
Radiant.set_dimension(this::Radiant.Geometry,dimension::Int64)
Radiant.set_boundary_conditions(this::Radiant.Geometry,boundary::String,boundary_condition::String)
Radiant.set_number_of_regions(this::Radiant.Geometry,axis::String,number_of_regions::Int64)
Radiant.set_voxels_per_region(this::Radiant.Geometry,axis::String,voxels_per_region::Vector{Int64})
Radiant.set_region_boundaries(this::Radiant.Geometry,axis::String,region_boundaries::Vector{Float64})

Radiant.set_material_per_region(this::Radiant.Geometry,material_per_region::Array{Radiant.Material})
Radiant.build(this::Radiant.Geometry,cs::Radiant.Cross_Sections)
```