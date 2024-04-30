# Transport calculations

RADIANT can compute the coupled transport of particles.

## Example 1: Coupled transport of photons-electrons-positrons for an electron beam normally incident on 1D water slab

The 'water' and 'cs' types come from [Example 1: Production of coupled photons-electrons-positrons cross-sections in water](example_cross_sections.md)

```julia
# Define geometry
geo = Geometry()
geo.set_type("Cartesian")
geo.set_dimension(1)
geo.set_material_per_region([water])
geo.set_boundary_conditions("X-","void")
geo.set_boundary_conditions("X+","void")
geo.set_number_of_regions("X",1)
geo.set_voxels_per_region("X",[80])
geo.set_region_boundaries("X",[0.0,10.0])
geo.build(cs)

# Define discretization methods for transport solver
m1 = Discrete_Ordinates()
m1.set_particle("photons")
m1.set_solver_type("BFP")
m1.set_acceleration("livolant")
m1.set_quadrature("Gauss-Lobatto",8)
m1.set_legendre_order(7)
m1.set_angular_boltzmann("galerkin-d")
m1.set_convergence_criterion(1e-5)
m1.set_maximum_iteration(200)
m1.set_scheme("x","DG",2)

m2 = Discrete_Ordinates()
m2.set_particle("electrons")
m2.set_solver_type("BFP")
m2.set_acceleration("livolant")
m2.set_quadrature("Gauss-Lobatto",8)
m2.set_legendre_order(7)
m2.set_angular_boltzmann("galerkin-d")
m2.set_angular_fokker_planck("finite-difference")
m2.set_convergence_criterion(1e-5)
m2.set_maximum_iteration(200)
m2.set_scheme("E","DG",2)
m2.set_scheme("x","DG",2)

m3 = Discrete_Ordinates()
m3.set_particle("positrons")
m3.set_solver_type("BFP")
m3.set_acceleration("livolant")
m3.set_quadrature("Gauss-Lobatto",8)
m3.set_legendre_order(7)
m3.set_angular_boltzmann("galerkin-d")
m3.set_angular_fokker_planck("finite-difference")
m3.set_convergence_criterion(1e-5)
m3.set_maximum_iteration(200)
m3.set_scheme("E","DG",2)
m3.set_scheme("x","DG",2)

solvers = Solvers()
solvers.add_solver(m1)
solvers.add_solver(m2)
solvers.add_solver(m3)
solvers.set_number_of_generations(2)

# Define fixed external sources
ss = Surface_Source()
ss.set_particle("electrons")
ss.set_intensity(1.0)
ss.set_energy_group(1)
ss.set_direction([1.0,0.0,0.0])
ss.set_location("X-")

s = Fixed_Sources(cs,geo,solvers)
s.add_source(ss)

# Consolidate the computation parameters and solve the transport problem
c1 = Computation_Unit()
c1.set_cross_sections(cs)
c1.set_geometry(geo)
c1.set_solvers(solvers)
c1.set_sources(s)
c1.run()

# Dose and energy spectrum solutions
x = c1.get_voxels_position("x")
E_γ = c1.get_energies("photons")
E_e = c1.get_energies("electrons")
E_p = c1.get_energies("positrons")
dose = c1.get_energy_deposition("total")
flux_γ = c1.get_flux("photons")
flux_e = c1.get_flux("electrons")
flux_p = c1.get_flux("positrons")

```