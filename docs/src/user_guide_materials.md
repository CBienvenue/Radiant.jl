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

!!! note "Additivity has a floor that `I` cannot lift"
    Any additivity-based stopping-power kernel misses NIST's liquid water by about
    28% at 1 keV and 18% at 100 keV. Applying additivity to NIST's *own* hydrogen and
    oxygen tables reproduces the same gap, so this is the assumption failing, not the
    implementation — ICRU 49 uses measured molecular data for water for exactly this
    reason. `set_mean_excitation_energy` does not fix it, because what fails is the
    stopping-power mixing rather than `I`. Only per-compound data would, which is what
    `set_stopping_power` (§2.6) is for.

## 2.6 Supplying a Stopping-Power Table

`set_stopping_power` reads the electronic stopping power of a material, for one
particle, out of a two-column text file instead of computing it. On this path neither
the corrected Bethe formula nor the SRIM tables are used: the file *is* the stopping
power.

This is the escape hatch, for a measured table, a table from a code Radiant does not
implement, or one a group has agreed on. It is attached to a `Material` rather than to
`Inelastic_Collision` because the table describes one material and nothing about it
generalises to the next one in the geometry.

```julia
water.set_stopping_power(Proton,"pstar_water_proton.txt")
water.get_stopping_power(Proton)   # resolved absolute path, or missing
water.get_stopping_power()         # the whole Dict
```

### File format

Plain text, two columns, conventionally placed in `data/`. A relative name is resolved
against that directory; an absolute path is used as given.

```
# proton in liquid water, PSTAR, retrieved 2026-01-14
# E [MeV]      S/rho [MeV cm^2/g]
1.000e-03      1.769e+02
1.500e-03      1.984e+02
...
1.000e+04      1.848e+00
```

Column 1 is the projectile's **total kinetic energy in MeV**, not energy per nucleon.
Column 2 is the **mass** stopping power in **MeV cm²/g**. The separator is whitespace or
a comma, blank lines are skipped, and a line starting with `#` or `%` is a comment.
Energies must be strictly increasing and are never sorted for you: an out-of-order row
is far more often a corrupt paste than a descending table.

A third column is an **error**, not something ignored. In SRIM's output that column is
*nuclear* stopping, and adding it in silently would be wrong by percent at the low end.

The file is parsed and validated inside `set_stopping_power`, so a malformed table is
reported at your own call site rather than deep inside `generate_cross_sections`.

### The table must be a total electronic stopping power

Radiant computes the soft stopping power as `Sr = Stot - Sc`, subtracting the
catastrophic part it derives from its own differential cross-section. So the table must
be:

  * **total**, not restricted and not soft-only — a restricted table subtracts the
    catastrophic part twice;
  * **electronic only** — a nuclear column, or the sum of the two, is the most likely
    thing to be pasted in by mistake, and below about 100 keV/u the nuclear term is
    percent-level and rising, small enough to survive a sanity check.

The soft/catastrophic split is therefore approximate on this path, in the same way it is
for the `"icru49"` low-energy branch: `Stot` comes from your source and `Sc` from
Radiant's analytic cross-section, and the two were never reconciled.

### One table, one particle

A table registered for `Proton` changes nothing for any other particle — electrons in
the same material still go through the computed path, and so do alphas. Radiant runs
coupled multi-particle problems, and a single un-keyed table would quietly apply a proton
table to electrons the moment an electron group was added. For two particles, call it
twice.

### What it overrides

For that one particle in that one material it supersedes `set_use_srim_sp`, `set_ecut`,
`set_emax`, `set_low_energy_model`, `set_low_energy_source`, `set_is_bragg_additivity`,
every correction inside the Bethe formula, and `set_mean_excitation_energy`. Those all
stay set and simply stop having an effect; none of the combinations is an error, because
handing in a table is switching sources rather than making a mistake.

It does **not** replace the radiative or nuclear stopping power. `Bremsstrahlung` and
`Elastic_Scattering` still contribute, so a Radiant total will differ from the table at
high energy.

### Density

The table is scaled by **this material's** density. For everything except the density
effect that is correct, since mass stopping power is density-independent — but the
density effect is baked into whatever numbers you supplied, at whatever density your
source assumed, and rescaling cannot touch it. A water-vapour table handed to a
liquid-water material gives vapour physics scaled to liquid density, silently. Radiant
cannot detect this and does not pretend to.

### Range

The table must span the whole energy range of the group structure for that particle. An
energy outside it raises an error naming the energy, the file and its range; it is never
extrapolated, clamped, or quietly replaced by the computed value, because a silent
fallback would put a discontinuity in the middle of a stopping-power curve where you
could not see it. The fix is always to extend the table or to narrow the group structure.

!!! warning "This is the one path with no validation behind it"
    The rest of Radiant's stopping-power chain was measured against NIST STAR and
    SRIM-2013 over hundreds of thousands of points. A user table is measured against
    nothing, because its correctness is yours to establish. Radiant validates the
    *file* — parsable, positive, strictly increasing, three rows or more — and then does
    exactly what it says. It does not compare your numbers against its own, and will not
    warn you when they disagree.

## 2.7 Inspecting a Material

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

## 2.8 Summary of the Material API

| Method                                   | Description                                                |
|------------------------------------------|------------------------------------------------------------|
| `Material(tag::String="")`               | Constructor with optional tag.                             |
| `set_tag(tag)`                           | Set the material tag.                                      |
| `set_density(ρ)`                         | Set the mass density [g/cm³].                              |
| `set_state_of_matter(s)`                 | Set the state of matter (`"solid"`, `"liquid"`, `"gaz"`).  |
| `add_element(sym,w=1)`                   | Add an element with its mass weight fraction.              |
| `set_mean_excitation_energy(I)`          | Override the mean excitation energy [eV] (else auto).      |
| `set_stopping_power(particle,file)`      | Read the electronic stopping power for one particle from a file. |
| `get_tag()`                              | Get the material tag.                                      |
| `get_density()`                          | Get the mass density.                                      |
| `get_state_of_matter()`                  | Get the state of matter.                                   |
| `get_number_of_elements()`               | Get the number of elements composing the material.         |
| `get_atomic_numbers()`                   | Get the atomic numbers of the elements.                    |
| `get_weight_fractions()`                 | Get the mass weight fractions of the elements.             |
| `get_mean_excitation_energy()`           | Get the user-defined mean excitation energy [eV] or `missing`. |
| `get_stopping_power(particle)`           | Get the registered stopping-power file for a particle, or `missing`. |
| `get_stopping_power()`                   | Get every registered stopping-power file.                  |
