# 2 Materials

Materials in Radiant are defined by a `Radiant.Material` object. A material is characterized by an optional human-readable tag, a mass density, a state of matter and a list of elements with their corresponding mass weight fractions.

## 2.1 Instantiating a Material

A new material is instantiated with the `Material` constructor, which accepts an optional string tag used to identify the material. This tag must be unique within a calculation if more than one material is defined.

```julia
water = Material("water")          # Tagged constructor (recommended)
unnamed = Material()               # Tag can also be set later with set_tag(...)
unnamed.set_tag("aluminum")
```

## 2.2 Predefined Materials

For convenience, Radiant ships with a large library of predefined materials based on the [NIST database](https://physics.nist.gov/PhysRefData/Star/Text/download.html) of elements, compounds and mixtures (water, air, lead, compact bone, A-150 tissue-equivalent plastic, and many more). Each one is a zero-argument constructor that returns a fully-configured `Material`, with its density, mean excitation energy, state of matter and elemental composition already set to the NIST-recommended values:

```julia
w   = Water()                       # liquid water
air = Air_Dry_Near_Sea_Level()      # dry air near sea level
pb  = Lead()                        # elemental lead
```

These are simply shortcuts for a material assembled by hand as in the next sections; any property can still be overridden afterwards (e.g. `w.set_density(0.998)`). The complete list of available constructors, together with the exact density, mean excitation energy and composition preset for each one, is given in the [Predefined Materials](@ref) page of the API reference.

## 2.3 Defining the Composition

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

## 2.4 State of Matter

By default a material is `"solid"`. The state of matter influences density-effect corrections in the inelastic-collision cross-sections. The accepted values are `"solid"`, `"liquid"` and `"gaz"`:

```julia
water.set_state_of_matter("liquid")
```

## 2.5 Mean Excitation Energy

The mean excitation energy ``I`` is the key material parameter of the collisional (inelastic) stopping power. By default Radiant determines it automatically:

1. for a few recognized compositions, the ICRU-recommended value is used (e.g. ``78`` eV for liquid water, ``85.7`` eV for dry air);
2. for any other material, ``I`` is computed from the elemental values through the Bragg additivity rule.

The additivity rule **systematically overestimates** ``I`` for many compounds because it ignores chemical-binding and phase effects (for example it yields ``≈ 87`` eV for an ICRP soft-tissue composition and ``≈ 92`` eV for striated muscle, versus the ICRU/ESTAR values of about ``73`` eV and ``74.7`` eV). Since ``I`` directly shifts the stopping power, this can move the Bragg peak noticeably. For tissues and other compounds it is therefore recommended to impose the ICRU/ESTAR value explicitly with `set_mean_excitation_energy` (in eV):

```julia
muscle.set_mean_excitation_energy(74.7)   # ICRU/ESTAR value for striated muscle
```

A value set this way takes precedence over both the tabulated value and the additivity rule. The getter `get_mean_excitation_energy()` returns the user-defined value, or `missing` if none was set (in which case the automatic value is used internally).

## 2.6 Inspecting a Material

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
water.get_mean_excitation_energy() # missing (78 eV used internally for water)
```

## 2.7 Summary of the Material API

| Method                                   | Description                                                |
|------------------------------------------|------------------------------------------------------------|
| `Material(tag::String="")`               | Constructor with optional tag.                             |
| `set_tag(tag)`                           | Set the material tag.                                      |
| `set_density(ρ)`                         | Set the mass density [g/cm³].                              |
| `set_state_of_matter(s)`                 | Set the state of matter (`"solid"`, `"liquid"`, `"gaz"`).  |
| `add_element(sym,w=1)`                   | Add an element with its mass weight fraction.              |
| `set_mean_excitation_energy(I)`          | Override the mean excitation energy [eV] (else auto).      |
| `get_tag()`                              | Get the material tag.                                      |
| `get_density()`                          | Get the mass density.                                      |
| `get_state_of_matter()`                  | Get the state of matter.                                   |
| `get_number_of_elements()`               | Get the number of elements composing the material.         |
| `get_atomic_numbers()`                   | Get the atomic numbers of the elements.                    |
| `get_weight_fractions()`                 | Get the mass weight fractions of the elements.             |
| `get_mean_excitation_energy()`           | Get the user-defined mean excitation energy [eV] or `missing`. |
