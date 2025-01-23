# 5 Geometry

A geometry in Radiant is represented by the `Radiant.Geometry` object. Once instantiated, its properties can be configured, and calculation can be done to prepare for transport calculations.

# 5.1 Cartesian geometries

Radiant provides the capabilities to produce heterogenous medium in Cartesian geometries. For 1D Cartesian geometry, only the axis `x` is defined, for 2D Cartesian geometry, only the axis `x` and `y`, and in 3D Cartesian geometry, the axis `x`, `y` and `z` are to be defined. Let's examine how to build a geometry object. The first step is to instantiate the `Radiant.Geometry` object, set the type of the geometry (Cartesian) and its dimension:

```julia
geo = Geometry()
geo.set_type("cartesian")
geo.set_dimension(1)
```

Then, the geometry has to be divided in regions, such as each region can be associated with a material. For exemle, in 1D geometry, a Al-Au-Al slab can be described as

```julia
geo.set_number_of_regions("x",3)                  # Set the number of regions along x
geo.set_region_boundaries("x",[0.0,0.2,0.4,0.6])  # Set the region boundaries.
geo.set_voxels_per_region("x",[20,20,20])         # Set the number of voxels per regions
geo.set_boundary_conditions("x-","void")          # Set the x- boundary condition
geo.set_boundary_conditions("x+","void")          # Set the x+ boundary condition
geo.set_material_per_region([al,au,al])           # Set the material in each region
```

Finally, the geometry component can be prepared for the further transport calculations using

```julia
geo.build(cs)
```

## 5.2 Other Kind of Geometries

!!! note
    At the moment, Radiant only treat 1D, 2D and 3D Cartesian geometries.