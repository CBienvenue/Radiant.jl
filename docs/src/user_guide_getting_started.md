# 1 Getting Started with Radiant

## 1.1 What is Radiant

Radiant is an open-source, general-purpose package for simulating the transport of ionizing radiation in matter using deterministic methods, in opposition to Monte-Carlo codes such as GEANT4 or EGSnrc. It is written in Julia, an open-source programming language providing Python-like readability and flexibility, combined with execution times comparable to those of C++ and FORTRAN [bezanson2017julia](@cite).

Radiant specializes in the coupled transport of photons, electrons and positrons with kinetic energies in the range 1 keV to 900 MeV, but also provides building blocks for heavy charged particles (protons, antiprotons, alphas, muons, antimuons) and user-defined custom particles.

## 1.2 How does it work

Radiant performs a calculation by first translating the user input into a set of data structures (also called *objects*) and then executing the transport solver on them. The user does not need to know the internal structure of the package: a calculation is built up by instantiating high-level objects and configuring their properties through *setter* methods (named `set_*`) and *action* methods (such as `build()`, `run()` or `add_*`). Once the calculation has run, results are retrieved through *getter* methods (named `get_*`).

The objects manipulated by the user are introduced one by one in the next sections of this User's Guide. They fall into the following categories:

| Object                  | Purpose                                                                       |
|-------------------------|-------------------------------------------------------------------------------|
| `Material`              | Defines a material (composition, density, state of matter).                   |
| `Particle`              | Represents a particle type (`Photon`, `Electron`, `Positron`, ...).           |
| `Interaction`           | Subtypes (`Compton`, `Bremsstrahlung`, ...) describe the physical processes.  |
| `Cross_Sections`        | Multigroup cross-section library, built or loaded from a file.                |
| `Geometry`              | Spatial layout (Cartesian 1D/2D/3D) and material map.                         |
| `SN`, `DPN`, `GN`       | Discretization methods, one per particle.                                     |
| `Solvers`               | Collection of discretization methods, controlling coupled transport.          |
| `Surface_Source`        | Directional source on a boundary of the geometry.                             |
| `Volume_Source`         | Isotropic source distributed inside the geometry.                             |
| `Fixed_Sources`         | Collection of surface and volume sources for a calculation.                   |
| `Computation_Unit`      | Assembles all of the above and executes the calculation.                      |

## 1.3 Method-call Notation

Throughout this guide, methods are called using the Python-like dot notation, e.g.

```julia
water.set_density(1.0)
```

This is equivalent to the functional notation

```julia
set_density(water,1.0)
```

Both forms work. The dot notation tends to read more naturally when an object is configured through several successive calls.

## 1.4 A First Complete Example

The following script gives the overall shape of a Radiant calculation. Each block corresponds to a section of this User's Guide.

```julia
using Radiant

# 1. Define materials
water = Material("water")
water.set_density(1.0)
water.set_state_of_matter("liquid")
water.add_element("H",0.1111)
water.add_element("O",0.8889)

# 2. Define particles
electron = Electron()
photon   = Photon()

# 3. Build the multigroup cross-sections library
cs = Cross_Sections()
cs.set_source("physics-models")
cs.set_materials([water])
cs.set_particles([electron,photon])
cs.set_group_structure("log",20,10.0,0.001)
cs.set_interactions([Inelastic_Collision(),Elastic_Collision(),Bremsstrahlung(),
                     Compton(),Photoelectric(),Pair_Production(),Annihilation(),
                     Rayleigh(),Relaxation()])
cs.set_legendre_order(7)
cs.build()

# 4. Define the geometry
geo = Geometry()
geo.set_type("cartesian")
geo.set_dimension(1)
geo.set_number_of_regions("x",1)
geo.set_region_boundaries("x",[0.0,5.0])
geo.set_voxels_per_region("x",[50])
geo.set_boundary_conditions("x-","void")
geo.set_boundary_conditions("x+","void")
geo.set_material_per_region([water])
geo.build(cs)

# 5. Define a discretization method per particle
sn_e = SN()
sn_e.set_particle(electron)
sn_e.set_solver_type("BFP")
sn_e.set_quadrature("gauss-lobatto",8)
sn_e.set_legendre_order(7)
sn_e.set_angular_boltzmann("galerkin")
sn_e.set_angular_fokker_planck("finite-difference")
sn_e.set_acceleration("livolant")
sn_e.set_scheme("x","DG",2)
sn_e.set_scheme("E","DG",2)

sn_p = SN()
sn_p.set_particle(photon)
sn_p.set_solver_type("BTE")
sn_p.set_quadrature("gauss-lobatto",8)
sn_p.set_legendre_order(7)
sn_p.set_scheme("x","DG",2)

solvers = Solvers()
solvers.add_solver(sn_e)
solvers.add_solver(sn_p)
solvers.set_number_of_generations(5)

# 6. Define the fixed sources
ss = Surface_Source()
ss.set_particle(photon)
ss.set_intensity(1.0)
ss.set_energy_group(1)
ss.set_direction([1.0,0.0,0.0])
ss.set_location("x-")

s = Fixed_Sources(cs,geo,solvers)
s.add_source(ss)

# 7. Run the calculation and post-process results
c = Computation_Unit()
c.set_cross_sections(cs)
c.set_geometry(geo)
c.set_solvers(solvers)
c.set_sources(s)
c.run()

x  = c.get_voxels_position("x")
de = c.get_energy_deposition()
dq = c.get_charge_deposition()
```

The remaining chapters of this User's Guide describe each step in detail, including the optional settings.

## 1.5 Running a Script as a Radiant Input File

A Julia script that uses Radiant can be executed directly from the command line. The `@radiant_input` macro can be added at the top of the script to enable this behavior:

```julia
using Radiant
@radiant_input

# ... script content ...
```

Such a file can then be run by passing its path to a Julia executable, or via `Radiant.run_script("path/to/script.jl")` from another Julia session.
