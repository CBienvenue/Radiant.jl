# Geometry

## Structure
```@docs
Radiant.Geometry
```

## Build
```@docs
Radiant.build(this::Radiant.Geometry,cs::Radiant.Cross_Sections)
```

## Setters
```@docs
Radiant.set_type(this::Radiant.Geometry,type::String)
Radiant.set_dimension(this::Radiant.Geometry,dimension::Int64)
Radiant.set_boundary_conditions(this::Radiant.Geometry,boundary::String,boundary_condition::String)
Radiant.set_number_of_regions(this::Radiant.Geometry,axis::String,number_of_regions::Int64)
Radiant.set_voxels_per_region(this::Radiant.Geometry,axis::String,voxels_per_region::Vector{Int64})
Radiant.set_region_boundaries(this::Radiant.Geometry,axis::String,region_boundaries::Vector{Float64})
Radiant.set_material_per_region(this::Radiant.Geometry,material_per_region::Array{Radiant.Material})
```

## Getters
```@docs
Radiant.get_type(this::Radiant.Geometry)
Radiant.get_dimension(this::Radiant.Geometry)
Radiant.get_axis(this::Radiant.Geometry)
Radiant.get_number_of_voxels(this::Radiant.Geometry)
Radiant.get_voxels_width(this::Radiant.Geometry)
Radiant.get_voxels_position(this::Radiant.Geometry,axis::String)
Radiant.get_voxels_position(this::Radiant.Geometry)
Radiant.get_voxels_boundaries(this::Radiant.Geometry,axis::String)
Radiant.get_voxels_boundaries(this::Radiant.Geometry)
Radiant.get_material_per_voxel(this::Radiant.Geometry)
Radiant.get_boundary_conditions(this::Radiant.Geometry,boundary::String)
Radiant.get_boundary_conditions(this::Radiant.Geometry)
```
