# 8 Electromagnetic Fields

Radiant can transport charged particles in the presence of an external, spatially-constant **magnetic field**, accounting for the Lorentz force that curves the particle trajectories. The field is represented by an `Electromagnetic_Field` object that is attached to the `Computation_Unit` alongside the cross-sections, geometry, solvers and sources.

## 8.1 The `Electromagnetic_Field` Object

An `Electromagnetic_Field` is instantiated and given a magnetic field vector, expressed in **tesla** along the `x`-, `y`- and `z`-axes of the geometry:

```julia
emf = Electromagnetic_Field()
emf.set_magnetic_field([0.0, 0.0, 1.5])   # (Bx, By, Bz) in tesla
```

It is then assigned to the `Computation_Unit`:

```julia
c.set_electromagnetic_field(emf)
```

The field is **optional**: if no `Electromagnetic_Field` is assigned, the calculation proceeds without one, exactly as in the previous chapters.

## 8.2 Requirements and Scope

The magnetic-field transport is subject to a few requirements:

- **Magnetic field only.** The electric field is not yet supported; assigning a non-zero electric field raises an error.
- **`SN` solver only.** The feature is implemented for the discrete-ordinates (`SN`) solver. Requesting a field with the `GN` solver raises an error.
- **A continuous-slowing-down solver type.** The field acts on the momentum of charged particles, so the solver type must be one of `"BFP"`, `"BCSD"`, `"CSD"` or `"FP"` (see Section 6).
- **A full-sphere angular quadrature.** Because the field couples every direction on the unit sphere, the quadrature must be defined over the whole sphere: `"lebedev"`, `"carlson"` or `"gauss-legendre-chebychev"`. This holds even in 1D and 2D geometries — the 1D-only `"gauss-legendre"` and `"gauss-lobatto"` quadratures cannot be used together with a field.

The field acts on each particle according to its electric charge; neutral particles (photons) are therefore unaffected. Kinetic energies are expressed in MeV and lengths in cm, consistently with the rest of Radiant, so that a magnetic field in tesla yields the physically correct gyration.

## 8.3 A Complete Example

The following script transports a 10 MeV electron beam into a water slab under a transverse magnetic field. The magnetic field bends the electron paths, which shifts the depth–dose curve relative to the field-free case.

```julia
using Radiant

# Material, particle and cross-sections
water = Water()
electron = Electron()
cs = Cross_Sections()
cs.set_source("physics-models")
cs.set_materials([water])
cs.set_particles([electron])
cs.set_group_structure("log",80,10.0,0.001)
cs.set_interactions([Elastic_Collision(),Inelastic_Collision(),Bremsstrahlung()])
cs.set_legendre_order(15)
cs.build()

# Geometry (1D water slab)
geo = Geometry()
geo.set_type("cartesian")
geo.set_dimension(1)
geo.set_number_of_regions("x",1)
geo.set_region_boundaries("x",[0.0,10.0])
geo.set_voxels_per_region("x",[80])
geo.set_boundary_conditions("x-","void")
geo.set_boundary_conditions("x+","void")
geo.set_material_per_region([water])
geo.build(cs)

# SN solver (BFP, full-sphere Lebedev quadrature)
sn = SN()
sn.set_particle(electron)
sn.set_solver_type("BFP")
sn.set_quadrature("lebedev",15)
sn.set_legendre_order(7)
sn.set_angular_boltzmann("galerkin")
sn.set_angular_fokker_planck("finite-difference")
sn.set_acceleration("gmres")
sn.set_scheme("x","DG",2)
sn.set_scheme("E","DG",2)
solvers = Solvers()
solvers.add_solver(sn)

# Electron beam entering along +x
ss = Surface_Source()
ss.set_particle(electron)
ss.set_intensity(1.0)
ss.set_energy_group(1)
ss.set_direction([1.0,0.0,0.0])
ss.set_location("x-")
s = Fixed_Sources(cs,geo,solvers)
s.add_source(ss)

# External magnetic field, transverse to the beam
emf = Electromagnetic_Field()
emf.set_magnetic_field([0.0,0.0,1.5])   # 1.5 T along z

# Assemble and run
c = Computation_Unit()
c.set_cross_sections(cs)
c.set_geometry(geo)
c.set_solvers(solvers)
c.set_sources(s)
c.set_electromagnetic_field(emf)
c.run()

x  = c.get_voxels_position("x")
de = c.get_energy_deposition()
```

Running the same calculation with and without the `c.set_electromagnetic_field(emf)` call produces the field and field-free depth–dose curves, respectively.

## 8.4 Summary of the `Electromagnetic_Field` API

| Method                              | Description                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `Electromagnetic_Field()`           | Constructor.                                                      |
| `set_magnetic_field([Bx,By,Bz])`    | Set the magnetic field vector [tesla] along the geometry axes.    |
| `get_magnetic_field()`              | Get the magnetic field vector.                                    |

The field is attached to a calculation through the `Computation_Unit`:

| Method                              | Description                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `set_electromagnetic_field(emf)`    | Assign the electromagnetic field to the `Computation_Unit`.       |
