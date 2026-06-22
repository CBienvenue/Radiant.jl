# Predefined Materials

Radiant ships with a library of ready-to-use materials based on the
[NIST database](https://physics.nist.gov/PhysRefData/Star/Text/download.html) of
elements, compounds and mixtures. Each entry below is a zero-argument constructor
that returns a fully-configured [`Radiant.Material`](@ref) object, with its mass
density, mean excitation energy ``I``, state of matter (when relevant) and
elemental composition (mass weight fractions) preset to the NIST-recommended
values.

These constructors are convenience shortcuts: the returned material is identical
to one assembled by hand with `Material(...)`, `set_density`,
`set_mean_excitation_energy` and `add_element`. Any property can be overridden
after construction.

```julia
using Radiant

w = Water()                          # liquid water, I = 78 eV, H/O composition
w.set_density(0.998)                 # override the preset density if needed

pb = Lead()                          # elemental lead, I = 823 eV, ρ = 11.35 g/cm³
```

The full list of available constructors (with the exact preset density, mean
excitation energy and composition of each material) follows.

```@autodocs
Modules = [Radiant]
Pages   = ["Material_List.jl"]
Order   = [:function]
```
