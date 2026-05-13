# 2 Materials

Materials in Radiant are defined by a `Radiant.Material` object. A material is characterized by an optional human-readable tag, a mass density, a state of matter and a list of elements with their corresponding mass weight fractions.

## 2.1 Instantiating a Material

A new material is instantiated with the `Material` constructor, which accepts an optional string tag used to identify the material. This tag must be unique within a calculation if more than one material is defined.

```julia
water = Material("water")          # Tagged constructor (recommended)
unnamed = Material()               # Tag can also be set later with set_tag(...)
unnamed.set_tag("aluminum")
```

## 2.2 Defining the Composition

Once a material is instantiated, elements are added with `add_element`. For each element, the chemical symbol and its mass weight fraction (between 0 and 1) are provided. The weight fractions across all the elements of a material must sum to at most 1.

```julia
water = Material("water")
water.set_density(1.0)             # Density in g/cm³
water.add_element("H",0.1111)      # Hydrogen with mass weight of 11.11 %
water.add_element("O",0.8889)      # Oxygen with mass weight of 88.89 %
```

For a monoelemental material, the weight fraction can be omitted; it defaults to 1.0 and the density is automatically set to the tabulated density of that element at 20°C:

```julia
al = Material("aluminum")
al.add_element("Al")               # 100 % Al, density set to 2.70 g/cm³
```

If `set_density(...)` is called later, it overrides this default value.

## 2.3 State of Matter

By default a material is `"solid"`. The state of matter influences density-effect corrections in the inelastic-collision cross-sections. The accepted values are `"solid"`, `"liquid"` and `"gaz"`:

```julia
water.set_state_of_matter("liquid")
```

## 2.4 Inspecting a Material

Calling `println(material)` prints a summary of the material:

```
Material:
   Density (g/cm³):              1.0
   Elements in the compound:     ["h", "o"]
   Weight fractions:             [0.1111, 0.8889]
   State of matter:              liquid
```

A series of getter methods is also provided to retrieve individual fields:

```julia
water.get_tag()                    # "water"
water.get_density()                # 1.0
water.get_state_of_matter()        # "liquid"
water.get_number_of_elements()     # 2
water.get_atomic_numbers()         # [1, 8]
water.get_weight_fractions()       # [0.1111, 0.8889]
```

## 2.5 Summary of the Material API

| Method                                   | Description                                                |
|------------------------------------------|------------------------------------------------------------|
| `Material(tag::String="")`               | Constructor with optional tag.                             |
| `set_tag(tag)`                           | Set the material tag.                                      |
| `set_density(ρ)`                         | Set the mass density [g/cm³].                              |
| `set_state_of_matter(s)`                 | Set the state of matter (`"solid"`, `"liquid"`, `"gaz"`).  |
| `add_element(sym,w=1)`                   | Add an element with its mass weight fraction.              |
| `get_tag()`                              | Get the material tag.                                      |
| `get_density()`                          | Get the mass density.                                      |
| `get_state_of_matter()`                  | Get the state of matter.                                   |
| `get_number_of_elements()`               | Get the number of elements composing the material.         |
| `get_atomic_numbers()`                   | Get the atomic numbers of the elements.                    |
| `get_weight_fractions()`                 | Get the mass weight fractions of the elements.             |
