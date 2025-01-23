# 7 Fixed External Sources

A fixed external source in Radiant is represented by the `Radiant.Fixed_Sources` object. Once instantiated, the user can add two kinds of source to it using

```julia
s = Fixed_Sources(cs,geo,solvers)     # Instantiation using cross-sections, geometry and solver object
s.add_source(surface_source)          # Adding a surface source
s.add_source(volume_source)           # Adding a volume source
```

where `surface_source` and `volume_source` are respectively `Radiant.Surface_Source` and `Radiant.Volume_Source` objects defined in the following sections.

## 7.1 Surface Sources

A surface source is represented by the `Radiant.Surface_Source` object. It is defined on a boundary of the geometry. For Cartesian geometry, these boundaries are identified as `x-`, `x+`, `y-`, `y+`, `z-` and `z+`. It should be noted that the surface source is associated with the given particle ID. In 1D Cartesian geometry, the definition of the surface source can take this form

```julia
ss = Surface_Source()
ss.set_particle(photon)            # Set the particle
ss.set_intensity(1.0)              # Set the source intensity
ss.set_energy_group(1)             # Set the energy group
ss.set_direction([1.0,0.0,0.0])    # Set the direction using direction cosine
ss.set_location("X-")              # Set location of the source
```

In 2D or 3D geometry, an additionnal term has to be added to define the size of the source on the given boundary. For example, in 2D Cartesian geometry, taking the same example, the user should add

```julia
ss.set_boundaries("y",[2.0,3.0])   # Set the boundaries of the source along axis y
```

## 7.2 Volume Sources

A surface source is represented by the `Radiant.Volume_Source` object. It is defined within the geometry domain. For 1D Cartesian geometry, a volume source can be defined as

```julia
vs = Volume_Source()
vs.set_particle(photon)
vs.set_intensity(1.0)
vs.set_energy_group(1)
vs.set_boundaries("x",[0.0,2.0])
```

and similarly for 2D and 3D Cartesian geometry, adding boundaries along `y` and `z` axis.