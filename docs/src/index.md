# Radiant

Radiant is a package that performs deterministic transport of ionizing radiation in matter. It specializes in the coupled transport of photons, electrons, and positrons with kinetic energies between 1 keV and 900 MeV.

## Installation

For [Installation Instructions](https://cbienvenue.github.io/Radiant.jl/quick_start/), please refer to the User's Guide.

## Overview of Features 

### Multigroup Cross-Sections

- **Production of multigroup cross-sections for ionizing radiation**
  - **Particles** : photons | electrons | positrons
  - **Interactions** : Photoelectric effect | Compton scattering | Pair production | Rayleigh scattering | Elastic scattering of leptons | Impact ionization | Bremsstrahlung | Annihilation | Relaxation (fluorescence and Auger electrons)
  - **Energy** : 1 keV to 900 MeV
  - **Energy discretization**: Linear | Logarithmic
  - Compound materials
- Custom cross-sections
- Read and write cross-sections files

### Deterministic Solver

- **Solver** : Boltzmann transport equation | Boltzmann Fokker-Planck equation | Coupled transport of particles
- **Medium** : Cartesian geometry in 1D, 2D and 3D | Heterogeneous medium | Void boundary conditions
- **Fixed external sources** : Isotropic volume sources | Monodirectional surface sources
- **Angular discretization**
  - Discrete ordinates method
  - Galerkin quadrature method
  - **Angular Fokker-Planck discretization** : Moment-preserving | Finite-difference
  - **Quadratures** : Gauss-Legendre | Gauss-Lobatto | Gauss-Chebychev product | Carlson's level-symmetric | Lebedev
- **Spatial discretization** : Generalized diamond difference schemes | Generalized discontinuous Galerkin schemes | Adaptive weighted schemes
- **Energy discretization** : Multigroup method | Continuous slowing-down discretization (same as spatial)
- **Acceleration methods** : Livolant acceleration

## Examples

Examples for the Radiant package can be found in the GitHub repository [Radiant - Examples](https://github.com/CBienvenue/Radiant-Examples).