# 7 Fixed External Sources

External sources injected into a Radiant calculation are represented by a `Fixed_Sources` object, which is a collection of `Surface_Source` and `Volume_Source` objects.

## 7.1 The `Fixed_Sources` Object

A `Fixed_Sources` object is instantiated by passing it the cross-sections, geometry, and solvers it must be consistent with. Sources are then added one at a time with `add_source`:

```julia
s = Fixed_Sources(cs,geo,solvers)
s.add_source(surface_source)         # surface_source :: Surface_Source
s.add_source(volume_source)          # volume_source  :: Volume_Source
```

When a `Computation_Unit` is run, the `Fixed_Sources` object is built automatically: the contributions of every source are mapped onto the multigroup angular-flux representation expected by the solver, and added to the appropriate spatial moments of the source term.

After a calculation, a normalization factor representing the sum of all source intensities can be retrieved with `s.get_normalization_factor()`.

## 7.2 Surface Sources

A `Surface_Source` represents a directional source applied to a face of the geometry. In Cartesian geometries, the available boundaries are `"x-"`, `"x+"`, `"y-"`, `"y+"`, `"z-"` and `"z+"`.

### 7.2.1 1D Cartesian Geometry

```julia
ss = Surface_Source()
ss.set_particle(photon)              # Particle to be emitted
ss.set_intensity(1.0)                # Intensity [#/cm^(N-1)], N is the geometry dimension
ss.set_energy_group(1)               # Energy group index in which particles are emitted
ss.set_direction([1.0,0.0,0.0])      # Direction cosines (μ,η,ξ); must be a unit vector
ss.set_location("x-")                # Face on which the source is located
```

In 1D, the source covers the whole face and only the previous calls are required.

### 7.2.2 2D and 3D Cartesian Geometries

In multidimensional geometries, the extent of the source on the chosen face must be set with `set_boundaries(axis,boundaries)`. The two boundary values are given in ascending order, in centimetres:

```julia
ss.set_boundaries("y",[2.0,3.0])     # 2D: extent along y
ss.set_boundaries("z",[0.5,1.5])     # 3D: extent along z (in addition to the y extent)
```

### 7.2.3 Optional Legendre Truncation

By default the boundary angular flux is expanded to a high Legendre order (64). The truncation can be reduced if a coarser representation of the directional source is acceptable:

```julia
ss.set_legendre_order(7)
```

### 7.2.4 Direction Cosines

`set_direction(v)` expects a three-component unit vector `[μ,η,ξ]`. The intensity is interpreted as the integral of the angular flux over the source area, multiplied by the direction cosine of the source with respect to the inward normal.

## 7.3 Volume Sources

A `Volume_Source` represents an isotropic source distributed within an axis-aligned box inside the geometry.

### 7.3.1 Setting up a Volume Source

```julia
vs = Volume_Source()
vs.set_particle(photon)
vs.set_intensity(1.0)                # Intensity [#/cm^N]
vs.set_energy_group(1)
vs.set_boundaries("x",[0.0,2.0])     # Always required (1D, 2D, 3D)
vs.set_boundaries("y",[0.0,2.0])     # Required in 2D and 3D
vs.set_boundaries("z",[0.0,2.0])     # Required in 3D
```

The `boundaries` for each axis are an ascending two-element vector giving the lower and upper bound of the source box along that axis, in centimetres.

## 7.4 Mixing Multiple Sources

Several surface and volume sources can be added to the same `Fixed_Sources` object. For example:

```julia
s = Fixed_Sources(cs,geo,solvers)

# Photon pencil beam entering from x-
ss = Surface_Source()
ss.set_particle(photon)
ss.set_energy_group(1)
ss.set_direction([1.0,0.0,0.0])
ss.set_location("x-")
s.add_source(ss)

# Distributed electron source in a region
vs = Volume_Source()
vs.set_particle(electron)
vs.set_energy_group(5)
vs.set_boundaries("x",[2.0,3.0])
s.add_source(vs)
```

Different particles, different energy groups and different intensities can be combined freely in a single calculation.

## 7.5 Summary of the Source API

### `Fixed_Sources`

| Method                                          | Description                                                    |
|-------------------------------------------------|----------------------------------------------------------------|
| `Fixed_Sources(cs,geo,solvers)`                 | Constructor.                                                   |
| `add_source(source)`                            | Append a surface or volume source.                             |
| `build()`                                       | Build the collection (called automatically by `Computation_Unit.run`). |
| `get_source(p)`                                 | Get the per-particle aggregated source.                        |
| `get_particles()`                               | List the particles with at least one source.                   |
| `get_normalization_factor()`                    | Sum of source intensities.                                     |

### `Surface_Source`

| Method                                          | Description                                                    |
|-------------------------------------------------|----------------------------------------------------------------|
| `Surface_Source()`                              | Constructor.                                                   |
| `set_particle(p)`                               | Particle emitted by the source.                                |
| `set_intensity(I)`                              | Source intensity.                                              |
| `set_energy_group(g)`                           | Index of the energy group of emission.                         |
| `set_direction([μ,η,ξ])`                        | Direction cosines (unit vector).                               |
| `set_location(face)`                            | Face on which the source is located.                           |
| `set_boundaries(axis,bnds)`                     | Extent of the source on the face (2D / 3D).                    |
| `set_legendre_order(L)`                         | Truncation order of the boundary flux expansion.               |

### `Volume_Source`

| Method                                          | Description                                                    |
|-------------------------------------------------|----------------------------------------------------------------|
| `Volume_Source()`                               | Constructor.                                                   |
| `set_particle(p)`                               | Particle emitted by the source.                                |
| `set_intensity(I)`                              | Source intensity.                                              |
| `set_energy_group(g)`                           | Index of the energy group of emission.                         |
| `set_boundaries(axis,bnds)`                     | Extent of the source box along each axis.                      |
