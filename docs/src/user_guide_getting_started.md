# 1 Getting Started with Radiant

## 1.1 What is Radiant

Radiant is an open-source and general-purpose package to simulate the transport of particles in matter based on deterministic calculations, in opposition to Monte-Carlo codes such a GEANT4 or EGSnrc. It is written in Julia, an open-source programming language providing Python-like readability and flexibility, combined with execution times comparable to those of C++ and FORTRAN [bezanson2017julia](@cite).

## 1.2 How does it works

In order to proceed to calculations, Radiant requires input from the user about the particles, the materials, the cross-sections, the discretization methods, the geometry, the sources, etc. It provides the user data structure (also called objects), which are defined in the following sections of the User's Guide, and corresponding functions (also called methods) to defines their properties. The user then does not have to understand the internal programming structure of Radiant to make it works.

