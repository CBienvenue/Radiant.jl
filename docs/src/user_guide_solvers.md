# 6 Solvers

In Radiant, solvers and their methodologies are represented by the `Radiant.Solvers` object. Once instantiated, the user can add solvers using

```julia
solvers = Solvers()
solvers.add_solver(sn)
```

where `sn` is a `Radiant.Discrete_Ordinates` objects, defined in the following section. Nonetheless, the user can specify a maximum number of particle generations to produce in coupled particle transport with

```julia
solvers.set_number_of_generations(3)
```

## 6.1 Discrete Ordinates Solvers

A discrete ordinates solver is represented by a `Radiant.Discrete_Ordinates` object. Once instantiated, the user associate it with a particle and defined its properties. For example, for 1D transport of photons using the Boltzmann transport equation (BTE), the user could write:

```julia
sn = Discrete_Ordinates()
sn.set_particle(photon)                            # Set particle
sn.set_solver_type("BTE")                          # Set the solver to be the BTE
sn.set_acceleration("livolant")                    # Accelerate using Livolant acceleration
sn.set_quadrature("gauss-lobatto",8)               # Set a 8-points Gauss-Lobatto quadrature
sn.set_legendre_order(7)                           # Set a Legendre truncation order
sn.set_angular_boltzmann("galerkin")               # Apply a Galerkin quadrature method
sn.set_convergence_criterion(4e-7)                 # Set a convergence criterion 
sn.set_maximum_iteration(300)                      # Set a maximum number of iterations
sn.set_scheme("x","AWD",1)                         # Set the spatial difference scheme
```

Another example, for 3D transport of electrons using the Boltzmann Fokker-Planck (BFP) equation, the user could write:

```julia
sn = Discrete_Ordinates()
sn.set_particle(electron)                          # Set particle
sn.set_solver_type("BFP")                          # Set the solver to be the BFP
sn.set_acceleration("livolant")                    # Accelerate using Livolant acceleration
sn.set_quadrature("lebedev",8)                     # Set a Lebedev quadrature
sn.set_legendre_order(15)                          # Set a Legendre truncation order
sn.set_angular_boltzmann("galerkin")               # Apply a Galerkin quadrature method
sn.set_angular_fokker_planck("finite-difference")  # Set angular Fokker-Planck scheme
sn.set_convergence_criterion(4e-7)                 # Set a convergence criterion 
sn.set_maximum_iteration(300)                      # Set a maximum number of iterations
sn.set_scheme("E","DG",2)                          # Set the spatial difference scheme
sn.set_scheme("x","DD",1)                          # Set the spatial difference scheme
sn.set_scheme("y","DD",1)                          # Set the spatial difference scheme
sn.set_scheme("z","DD",1)                          # Set the spatial difference scheme
```

## 6.2 Other Kind of Solvers

!!! note
    At the moment, Radiant only treat discrete ordinates solvers.

