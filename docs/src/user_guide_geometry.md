# 5 Geometry

A geometry in Radiant is represented by the `Radiant.Geometry` object. Once instantiated, its dimension, region layout, voxelization, material map and boundary conditions are configured, and the structure is finalized with a call to `build(...)`.

## 5.1 Cartesian Geometries

Radiant currently treats Cartesian geometries in 1D, 2D and 3D. The geometry is described region by region along each axis:

- In 1D the only axis is `"x"`.
- In 2D the axes are `"x"` and `"y"`.
- In 3D the axes are `"x"`, `"y"` and `"z"`.

### 5.1.1 Instantiating and Setting the Dimension

```julia
geo = Geometry()
geo.set_type("cartesian")
geo.set_dimension(1)        # 1D Cartesian geometry
```

`set_dimension` only accepts 1, 2 or 3.

### 5.1.2 Defining Regions and Voxelization

For each axis, the user specifies:
1. The number of regions, with `set_number_of_regions(axis,N)`.
2. The region boundaries, in ascending order, with `set_region_boundaries(axis,bnds)`. The length of the boundary vector must be `N+1`.
3. The number of voxels used to discretize each region, with `set_voxels_per_region(axis,nvox)`. The length of `nvox` must be `N`.

For a 1D slab of three regions (Al, Au, Al) sharing the same voxel count:

```julia
geo.set_number_of_regions("x",3)
geo.set_region_boundaries("x",[0.0,0.2,0.4,0.6])
geo.set_voxels_per_region("x",[20,20,20])
```

The same pattern is used along `"y"` and `"z"` in higher dimensions.

### 5.1.3 Material Map

`set_material_per_region` takes a multi-dimensional array of `Material` objects, with one entry per region:

```julia
# 1D: a row vector of length Nx
geo.set_material_per_region([al,au,al])

# 2D: a matrix of size (Nx, Ny)
geo.set_material_per_region([al au al ; au al au])

# 3D: a 3-D array of size (Nx, Ny, Nz)
```

The materials placed in this array must be among those that were given to the corresponding `Cross_Sections` object.

### 5.1.4 Boundary Conditions

A boundary condition is set independently on each face of the domain. The face is identified by a string `"x-"`, `"x+"`, `"y-"`, `"y+"`, `"z-"`, `"z+"`, and the condition can be `"void"`, `"reflective"` or `"periodic"`:

```julia
geo.set_boundary_conditions("x-","void")
geo.set_boundary_conditions("x+","void")
```

In 2D and 3D, the conditions on the additional boundaries must also be set.

### 5.1.5 Building the Geometry

After the inputs are complete, the geometry is finalized with:

```julia
geo.build(cs)
```

`build(cs)` validates the inputs and pre-computes voxel widths, midpoint positions, boundaries, the material index in each voxel and the voxel volumes. The `Cross_Sections` argument is used to match the material objects from the geometry to those stored in the library.

## 5.2 Inspecting the Geometry

After `build(...)`, the following getters expose the discretization:

```julia
geo.get_type()                          # "cartesian"
geo.get_dimension()                     # 1, 2 or 3
geo.get_axis()                          # ["x"], ["x","y"], or ["x","y","z"]
geo.get_number_of_voxels()              # Vector with the number of voxels along each axis
geo.get_voxels_width()                  # Voxel widths along each axis
geo.get_voxels_position("x")            # Midpoint positions along the x-axis
geo.get_voxels_boundaries("x")          # Voxel boundaries along the x-axis
geo.get_material_per_voxel()            # Material index in each voxel
geo.get_boundary_conditions()           # All boundary conditions as integer codes
geo.get_boundary_conditions("x-")       # Boundary condition string at x-
```

## 5.3 Complete Example (1D Al-Au-Al slab)

```julia
geo = Geometry()
geo.set_type("cartesian")
geo.set_dimension(1)
geo.set_number_of_regions("x",3)
geo.set_region_boundaries("x",[0.0,0.2,0.4,0.6])
geo.set_voxels_per_region("x",[20,20,20])
geo.set_boundary_conditions("x-","void")
geo.set_boundary_conditions("x+","void")
geo.set_material_per_region([al,au,al])
geo.build(cs)
```

## 5.4 Other Kinds of Geometries

!!! note
    At the moment, Radiant only treats 1D, 2D and 3D Cartesian geometries. Curvilinear (cylindrical, spherical) and unstructured geometries are not yet supported.

## 5.5 Summary of the Geometry API

| Method                                          | Description                                                  |
|-------------------------------------------------|--------------------------------------------------------------|
| `Geometry()`                                    | Constructor.                                                 |
| `set_type(t)`                                   | Set the geometry type (`"cartesian"`).                       |
| `set_dimension(d)`                              | Set the geometry dimension (1, 2 or 3).                      |
| `set_number_of_regions(axis,N)`                 | Set the number of regions along an axis.                     |
| `set_region_boundaries(axis,bnds)`              | Set the region boundaries along an axis.                     |
| `set_voxels_per_region(axis,nvox)`              | Set the number of voxels per region along an axis.           |
| `set_material_per_region(mat_array)`            | Set the material in each region.                             |
| `set_boundary_conditions(face,bc)`              | Set the boundary condition on a given face.                  |
| `build(cs)`                                     | Build the discretized geometry.                              |
| `get_type()`, `get_dimension()`, `get_axis()`   | Geometry metadata.                                           |
| `get_number_of_voxels()`                        | Number of voxels along each axis.                            |
| `get_voxels_width()`                            | Voxel widths along each axis.                                |
| `get_voxels_position(axis)`                     | Midpoint positions along a given axis.                       |
| `get_voxels_boundaries(axis)`                   | Voxel boundaries along a given axis.                         |
| `get_material_per_voxel()`                      | Material index assigned to each voxel.                       |
| `get_boundary_conditions([face])`               | Get the boundary conditions (one face or all).               |
