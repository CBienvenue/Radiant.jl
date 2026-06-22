#----------------------------------------------------------------#
#
# This file defines a library of predefined-material constructors based on the
# NIST database of elements, compounds and mixtures. Each zero-argument function
# (e.g. `Hydrogen()`, `Water()`, `Bone_Compact_Icru()`) builds and returns a
# fully-configured `Material` object, with its mass density, mean excitation
# energy, state of matter (when relevant) and elemental composition (mass weight
# fractions) preset to the NIST-recommended values. They are convenience
# shortcuts: the returned `Material` is identical to one assembled by hand with
# `Material(...)`, `set_density`, `set_mean_excitation_energy` and `add_element`,
# and any property can be overridden afterwards.
#
# All constructors listed in `MATERIAL_CONSTRUCTORS` (at the bottom of this file)
# are exported by the Radiant module.
#
# References:
# - National Institute of Standards and Technology (NIST) database: https://physics.nist.gov/PhysRefData/Star/Text/download.html.
#
#----------------------------------------------------------------#

"""
    Hydrogen()

Build the predefined NIST material **HYDROGEN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.3748e-05` g/cm³
- Mean excitation energy ``I``: `19.2` eV
- Composition (element ⇒ mass weight fraction): `H` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Hydrogen()
    mat = Material("hydrogen")
    mat.set_density(8.3748e-05)
    mat.set_mean_excitation_energy(19.2) # eV
    mat.add_element("H",1.000000)
    return mat
end

"""
    Helium()

Build the predefined NIST material **HELIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.000166322` g/cm³
- Mean excitation energy ``I``: `41.8` eV
- Composition (element ⇒ mass weight fraction): `He` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Helium()
    mat = Material("helium")
    mat.set_density(0.000166322)
    mat.set_mean_excitation_energy(41.8) # eV
    mat.add_element("He",1.000000)
    return mat
end

"""
    Lithium()

Build the predefined NIST material **LITHIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.534` g/cm³
- Mean excitation energy ``I``: `40` eV
- Composition (element ⇒ mass weight fraction): `Li` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium()
    mat = Material("lithium")
    mat.set_density(0.534)
    mat.set_mean_excitation_energy(40) # eV
    mat.add_element("Li",1.000000)
    return mat
end

"""
    Beryllium()

Build the predefined NIST material **BERYLLIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.848` g/cm³
- Mean excitation energy ``I``: `63.7` eV
- Composition (element ⇒ mass weight fraction): `Be` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Beryllium()
    mat = Material("beryllium")
    mat.set_density(1.848)
    mat.set_mean_excitation_energy(63.7) # eV
    mat.add_element("Be",1.000000)
    return mat
end

"""
    Boron()

Build the predefined NIST material **BORON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.37` g/cm³
- Mean excitation energy ``I``: `76` eV
- Composition (element ⇒ mass weight fraction): `B` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Boron()
    mat = Material("boron")
    mat.set_density(2.37)
    mat.set_mean_excitation_energy(76) # eV
    mat.add_element("B",1.000000)
    return mat
end

"""
    Amorphous_Carbon()

Build the predefined NIST material **AMORPHOUS CARBON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2` g/cm³
- Mean excitation energy ``I``: `81` eV
- Composition (element ⇒ mass weight fraction): `C` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Amorphous_Carbon()
    mat = Material("amorphous_carbon")
    mat.set_density(2)
    mat.set_mean_excitation_energy(81) # eV
    mat.add_element("C",1.000000)
    return mat
end

"""
    Nitrogen()

Build the predefined NIST material **NITROGEN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00116528` g/cm³
- Mean excitation energy ``I``: `82` eV
- Composition (element ⇒ mass weight fraction): `N` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nitrogen()
    mat = Material("nitrogen")
    mat.set_density(0.00116528)
    mat.set_mean_excitation_energy(82) # eV
    mat.add_element("N",1.000000)
    return mat
end

"""
    Oxygen()

Build the predefined NIST material **OXYGEN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00133151` g/cm³
- Mean excitation energy ``I``: `95` eV
- Composition (element ⇒ mass weight fraction): `O` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Oxygen()
    mat = Material("oxygen")
    mat.set_density(0.00133151)
    mat.set_mean_excitation_energy(95) # eV
    mat.add_element("O",1.000000)
    return mat
end

"""
    Fluorine()

Build the predefined NIST material **FLUORINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00158029` g/cm³
- Mean excitation energy ``I``: `115` eV
- Composition (element ⇒ mass weight fraction): `F` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Fluorine()
    mat = Material("fluorine")
    mat.set_density(0.00158029)
    mat.set_mean_excitation_energy(115) # eV
    mat.add_element("F",1.000000)
    return mat
end

"""
    Neon()

Build the predefined NIST material **NEON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.000838505` g/cm³
- Mean excitation energy ``I``: `137` eV
- Composition (element ⇒ mass weight fraction): `Ne` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Neon()
    mat = Material("neon")
    mat.set_density(0.000838505)
    mat.set_mean_excitation_energy(137) # eV
    mat.add_element("Ne",1.000000)
    return mat
end

"""
    Sodium()

Build the predefined NIST material **SODIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.971` g/cm³
- Mean excitation energy ``I``: `149` eV
- Composition (element ⇒ mass weight fraction): `Na` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sodium()
    mat = Material("sodium")
    mat.set_density(0.971)
    mat.set_mean_excitation_energy(149) # eV
    mat.add_element("Na",1.000000)
    return mat
end

"""
    Magnesium()

Build the predefined NIST material **MAGNESIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.74` g/cm³
- Mean excitation energy ``I``: `156` eV
- Composition (element ⇒ mass weight fraction): `Mg` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Magnesium()
    mat = Material("magnesium")
    mat.set_density(1.74)
    mat.set_mean_excitation_energy(156) # eV
    mat.add_element("Mg",1.000000)
    return mat
end

"""
    Aluminum()

Build the predefined NIST material **ALUMINUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.6989` g/cm³
- Mean excitation energy ``I``: `166` eV
- Composition (element ⇒ mass weight fraction): `Al` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Aluminum()
    mat = Material("aluminum")
    mat.set_density(2.6989)
    mat.set_mean_excitation_energy(166) # eV
    mat.add_element("Al",1.000000)
    return mat
end

"""
    Silicon()

Build the predefined NIST material **SILICON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.33` g/cm³
- Mean excitation energy ``I``: `173` eV
- Composition (element ⇒ mass weight fraction): `Si` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silicon()
    mat = Material("silicon")
    mat.set_density(2.33)
    mat.set_mean_excitation_energy(173) # eV
    mat.add_element("Si",1.000000)
    return mat
end

"""
    Phosphorus()

Build the predefined NIST material **PHOSPHORUS** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.2` g/cm³
- Mean excitation energy ``I``: `173` eV
- Composition (element ⇒ mass weight fraction): `P` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Phosphorus()
    mat = Material("phosphorus")
    mat.set_density(2.2)
    mat.set_mean_excitation_energy(173) # eV
    mat.add_element("P",1.000000)
    return mat
end

"""
    Sulfur()

Build the predefined NIST material **SULFUR** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2` g/cm³
- Mean excitation energy ``I``: `180` eV
- Composition (element ⇒ mass weight fraction): `S` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sulfur()
    mat = Material("sulfur")
    mat.set_density(2)
    mat.set_mean_excitation_energy(180) # eV
    mat.add_element("S",1.000000)
    return mat
end

"""
    Chlorine()

Build the predefined NIST material **CHLORINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00299473` g/cm³
- Mean excitation energy ``I``: `174` eV
- Composition (element ⇒ mass weight fraction): `Cl` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Chlorine()
    mat = Material("chlorine")
    mat.set_density(0.00299473)
    mat.set_mean_excitation_energy(174) # eV
    mat.add_element("Cl",1.000000)
    return mat
end

"""
    Argon()

Build the predefined NIST material **ARGON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00166201` g/cm³
- Mean excitation energy ``I``: `188` eV
- Composition (element ⇒ mass weight fraction): `Ar` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Argon()
    mat = Material("argon")
    mat.set_density(0.00166201)
    mat.set_mean_excitation_energy(188) # eV
    mat.add_element("Ar",1.000000)
    return mat
end

"""
    Potassium()

Build the predefined NIST material **POTASSIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.862` g/cm³
- Mean excitation energy ``I``: `190` eV
- Composition (element ⇒ mass weight fraction): `K` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Potassium()
    mat = Material("potassium")
    mat.set_density(0.862)
    mat.set_mean_excitation_energy(190) # eV
    mat.add_element("K",1.000000)
    return mat
end

"""
    Calcium()

Build the predefined NIST material **CALCIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.55` g/cm³
- Mean excitation energy ``I``: `191` eV
- Composition (element ⇒ mass weight fraction): `Ca` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium()
    mat = Material("calcium")
    mat.set_density(1.55)
    mat.set_mean_excitation_energy(191) # eV
    mat.add_element("Ca",1.000000)
    return mat
end

"""
    Scandium()

Build the predefined NIST material **SCANDIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.989` g/cm³
- Mean excitation energy ``I``: `216` eV
- Composition (element ⇒ mass weight fraction): `Sc` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Scandium()
    mat = Material("scandium")
    mat.set_density(2.989)
    mat.set_mean_excitation_energy(216) # eV
    mat.add_element("Sc",1.000000)
    return mat
end

"""
    Titanium()

Build the predefined NIST material **TITANIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.54` g/cm³
- Mean excitation energy ``I``: `233` eV
- Composition (element ⇒ mass weight fraction): `Ti` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Titanium()
    mat = Material("titanium")
    mat.set_density(4.54)
    mat.set_mean_excitation_energy(233) # eV
    mat.add_element("Ti",1.000000)
    return mat
end

"""
    Vanadium()

Build the predefined NIST material **VANADIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.11` g/cm³
- Mean excitation energy ``I``: `245` eV
- Composition (element ⇒ mass weight fraction): `V` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Vanadium()
    mat = Material("vanadium")
    mat.set_density(6.11)
    mat.set_mean_excitation_energy(245) # eV
    mat.add_element("V",1.000000)
    return mat
end

"""
    Chromium()

Build the predefined NIST material **CHROMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.18` g/cm³
- Mean excitation energy ``I``: `257` eV
- Composition (element ⇒ mass weight fraction): `Cr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Chromium()
    mat = Material("chromium")
    mat.set_density(7.18)
    mat.set_mean_excitation_energy(257) # eV
    mat.add_element("Cr",1.000000)
    return mat
end

"""
    Manganese()

Build the predefined NIST material **MANGANESE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.44` g/cm³
- Mean excitation energy ``I``: `272` eV
- Composition (element ⇒ mass weight fraction): `Mn` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Manganese()
    mat = Material("manganese")
    mat.set_density(7.44)
    mat.set_mean_excitation_energy(272) # eV
    mat.add_element("Mn",1.000000)
    return mat
end

"""
    Iron()

Build the predefined NIST material **IRON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.874` g/cm³
- Mean excitation energy ``I``: `286` eV
- Composition (element ⇒ mass weight fraction): `Fe` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Iron()
    mat = Material("iron")
    mat.set_density(7.874)
    mat.set_mean_excitation_energy(286) # eV
    mat.add_element("Fe",1.000000)
    return mat
end

"""
    Cobalt()

Build the predefined NIST material **COBALT** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.9` g/cm³
- Mean excitation energy ``I``: `297` eV
- Composition (element ⇒ mass weight fraction): `Co` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cobalt()
    mat = Material("cobalt")
    mat.set_density(8.9)
    mat.set_mean_excitation_energy(297) # eV
    mat.add_element("Co",1.000000)
    return mat
end

"""
    Nickel()

Build the predefined NIST material **NICKEL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.902` g/cm³
- Mean excitation energy ``I``: `311` eV
- Composition (element ⇒ mass weight fraction): `Ni` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nickel()
    mat = Material("nickel")
    mat.set_density(8.902)
    mat.set_mean_excitation_energy(311) # eV
    mat.add_element("Ni",1.000000)
    return mat
end

"""
    Copper()

Build the predefined NIST material **COPPER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.96` g/cm³
- Mean excitation energy ``I``: `322` eV
- Composition (element ⇒ mass weight fraction): `Cu` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Copper()
    mat = Material("copper")
    mat.set_density(8.96)
    mat.set_mean_excitation_energy(322) # eV
    mat.add_element("Cu",1.000000)
    return mat
end

"""
    Zinc()

Build the predefined NIST material **ZINC** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.133` g/cm³
- Mean excitation energy ``I``: `330` eV
- Composition (element ⇒ mass weight fraction): `Zn` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Zinc()
    mat = Material("zinc")
    mat.set_density(7.133)
    mat.set_mean_excitation_energy(330) # eV
    mat.add_element("Zn",1.000000)
    return mat
end

"""
    Gallium()

Build the predefined NIST material **GALLIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.904` g/cm³
- Mean excitation energy ``I``: `334` eV
- Composition (element ⇒ mass weight fraction): `Ga` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gallium()
    mat = Material("gallium")
    mat.set_density(5.904)
    mat.set_mean_excitation_energy(334) # eV
    mat.add_element("Ga",1.000000)
    return mat
end

"""
    Germanium()

Build the predefined NIST material **GERMANIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.323` g/cm³
- Mean excitation energy ``I``: `350` eV
- Composition (element ⇒ mass weight fraction): `Ge` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Germanium()
    mat = Material("germanium")
    mat.set_density(5.323)
    mat.set_mean_excitation_energy(350) # eV
    mat.add_element("Ge",1.000000)
    return mat
end

"""
    Arsenic()

Build the predefined NIST material **ARSENIC** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.73` g/cm³
- Mean excitation energy ``I``: `347` eV
- Composition (element ⇒ mass weight fraction): `As` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Arsenic()
    mat = Material("arsenic")
    mat.set_density(5.73)
    mat.set_mean_excitation_energy(347) # eV
    mat.add_element("As",1.000000)
    return mat
end

"""
    Selenium()

Build the predefined NIST material **SELENIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.5` g/cm³
- Mean excitation energy ``I``: `348` eV
- Composition (element ⇒ mass weight fraction): `Se` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Selenium()
    mat = Material("selenium")
    mat.set_density(4.5)
    mat.set_mean_excitation_energy(348) # eV
    mat.add_element("Se",1.000000)
    return mat
end

"""
    Bromine()

Build the predefined NIST material **BROMINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00707218` g/cm³
- Mean excitation energy ``I``: `343` eV
- Composition (element ⇒ mass weight fraction): `Br` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bromine()
    mat = Material("bromine")
    mat.set_density(0.00707218)
    mat.set_mean_excitation_energy(343) # eV
    mat.add_element("Br",1.000000)
    return mat
end

"""
    Krypton()

Build the predefined NIST material **KRYPTON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00347832` g/cm³
- Mean excitation energy ``I``: `352` eV
- Composition (element ⇒ mass weight fraction): `Kr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Krypton()
    mat = Material("krypton")
    mat.set_density(0.00347832)
    mat.set_mean_excitation_energy(352) # eV
    mat.add_element("Kr",1.000000)
    return mat
end

"""
    Rubidium()

Build the predefined NIST material **RUBIDIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.532` g/cm³
- Mean excitation energy ``I``: `363` eV
- Composition (element ⇒ mass weight fraction): `Rb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rubidium()
    mat = Material("rubidium")
    mat.set_density(1.532)
    mat.set_mean_excitation_energy(363) # eV
    mat.add_element("Rb",1.000000)
    return mat
end

"""
    Strontium()

Build the predefined NIST material **STRONTIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.54` g/cm³
- Mean excitation energy ``I``: `366` eV
- Composition (element ⇒ mass weight fraction): `Sr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Strontium()
    mat = Material("strontium")
    mat.set_density(2.54)
    mat.set_mean_excitation_energy(366) # eV
    mat.add_element("Sr",1.000000)
    return mat
end

"""
    Yttrium()

Build the predefined NIST material **YTTRIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.469` g/cm³
- Mean excitation energy ``I``: `379` eV
- Composition (element ⇒ mass weight fraction): `Y` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Yttrium()
    mat = Material("yttrium")
    mat.set_density(4.469)
    mat.set_mean_excitation_energy(379) # eV
    mat.add_element("Y",1.000000)
    return mat
end

"""
    Zirconium()

Build the predefined NIST material **ZIRCONIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.506` g/cm³
- Mean excitation energy ``I``: `393` eV
- Composition (element ⇒ mass weight fraction): `Zr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Zirconium()
    mat = Material("zirconium")
    mat.set_density(6.506)
    mat.set_mean_excitation_energy(393) # eV
    mat.add_element("Zr",1.000000)
    return mat
end

"""
    Niobium()

Build the predefined NIST material **NIOBIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.57` g/cm³
- Mean excitation energy ``I``: `417` eV
- Composition (element ⇒ mass weight fraction): `Nb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Niobium()
    mat = Material("niobium")
    mat.set_density(8.57)
    mat.set_mean_excitation_energy(417) # eV
    mat.add_element("Nb",1.000000)
    return mat
end

"""
    Molybdenum()

Build the predefined NIST material **MOLYBDENUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `10.22` g/cm³
- Mean excitation energy ``I``: `424` eV
- Composition (element ⇒ mass weight fraction): `Mo` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Molybdenum()
    mat = Material("molybdenum")
    mat.set_density(10.22)
    mat.set_mean_excitation_energy(424) # eV
    mat.add_element("Mo",1.000000)
    return mat
end

"""
    Technetium()

Build the predefined NIST material **TECHNETIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.5` g/cm³
- Mean excitation energy ``I``: `428` eV
- Composition (element ⇒ mass weight fraction): `Tc` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Technetium()
    mat = Material("technetium")
    mat.set_density(11.5)
    mat.set_mean_excitation_energy(428) # eV
    mat.add_element("Tc",1.000000)
    return mat
end

"""
    Ruthenium()

Build the predefined NIST material **RUTHENIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `12.41` g/cm³
- Mean excitation energy ``I``: `441` eV
- Composition (element ⇒ mass weight fraction): `Ru` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ruthenium()
    mat = Material("ruthenium")
    mat.set_density(12.41)
    mat.set_mean_excitation_energy(441) # eV
    mat.add_element("Ru",1.000000)
    return mat
end

"""
    Rhodium()

Build the predefined NIST material **RHODIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `12.41` g/cm³
- Mean excitation energy ``I``: `449` eV
- Composition (element ⇒ mass weight fraction): `Rh` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rhodium()
    mat = Material("rhodium")
    mat.set_density(12.41)
    mat.set_mean_excitation_energy(449) # eV
    mat.add_element("Rh",1.000000)
    return mat
end

"""
    Palladium()

Build the predefined NIST material **PALLADIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `12.02` g/cm³
- Mean excitation energy ``I``: `470` eV
- Composition (element ⇒ mass weight fraction): `Pd` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Palladium()
    mat = Material("palladium")
    mat.set_density(12.02)
    mat.set_mean_excitation_energy(470) # eV
    mat.add_element("Pd",1.000000)
    return mat
end

"""
    Silver()

Build the predefined NIST material **SILVER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `10.5` g/cm³
- Mean excitation energy ``I``: `470` eV
- Composition (element ⇒ mass weight fraction): `Ag` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silver()
    mat = Material("silver")
    mat.set_density(10.5)
    mat.set_mean_excitation_energy(470) # eV
    mat.add_element("Ag",1.000000)
    return mat
end

"""
    Cadmium()

Build the predefined NIST material **CADMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.65` g/cm³
- Mean excitation energy ``I``: `469` eV
- Composition (element ⇒ mass weight fraction): `Cd` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cadmium()
    mat = Material("cadmium")
    mat.set_density(8.65)
    mat.set_mean_excitation_energy(469) # eV
    mat.add_element("Cd",1.000000)
    return mat
end

"""
    Indium()

Build the predefined NIST material **INDIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.31` g/cm³
- Mean excitation energy ``I``: `488` eV
- Composition (element ⇒ mass weight fraction): `In` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Indium()
    mat = Material("indium")
    mat.set_density(7.31)
    mat.set_mean_excitation_energy(488) # eV
    mat.add_element("In",1.000000)
    return mat
end

"""
    Tin()

Build the predefined NIST material **TIN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.31` g/cm³
- Mean excitation energy ``I``: `488` eV
- Composition (element ⇒ mass weight fraction): `Sn` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tin()
    mat = Material("tin")
    mat.set_density(7.31)
    mat.set_mean_excitation_energy(488) # eV
    mat.add_element("Sn",1.000000)
    return mat
end

"""
    Antimony()

Build the predefined NIST material **ANTIMONY** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.691` g/cm³
- Mean excitation energy ``I``: `487` eV
- Composition (element ⇒ mass weight fraction): `Sb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Antimony()
    mat = Material("antimony")
    mat.set_density(6.691)
    mat.set_mean_excitation_energy(487) # eV
    mat.add_element("Sb",1.000000)
    return mat
end

"""
    Tellurium()

Build the predefined NIST material **TELLURIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.24` g/cm³
- Mean excitation energy ``I``: `485` eV
- Composition (element ⇒ mass weight fraction): `Te` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tellurium()
    mat = Material("tellurium")
    mat.set_density(6.24)
    mat.set_mean_excitation_energy(485) # eV
    mat.add_element("Te",1.000000)
    return mat
end

"""
    Iodine()

Build the predefined NIST material **IODINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.93` g/cm³
- Mean excitation energy ``I``: `491` eV
- Composition (element ⇒ mass weight fraction): `I` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Iodine()
    mat = Material("iodine")
    mat.set_density(4.93)
    mat.set_mean_excitation_energy(491) # eV
    mat.add_element("I",1.000000)
    return mat
end

"""
    Xenon()

Build the predefined NIST material **XENON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00548536` g/cm³
- Mean excitation energy ``I``: `482` eV
- Composition (element ⇒ mass weight fraction): `Xe` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Xenon()
    mat = Material("xenon")
    mat.set_density(0.00548536)
    mat.set_mean_excitation_energy(482) # eV
    mat.add_element("Xe",1.000000)
    return mat
end

"""
    Cesium()

Build the predefined NIST material **CESIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.873` g/cm³
- Mean excitation energy ``I``: `488` eV
- Composition (element ⇒ mass weight fraction): `Cs` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cesium()
    mat = Material("cesium")
    mat.set_density(1.873)
    mat.set_mean_excitation_energy(488) # eV
    mat.add_element("Cs",1.000000)
    return mat
end

"""
    Barium()

Build the predefined NIST material **BARIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.5` g/cm³
- Mean excitation energy ``I``: `491` eV
- Composition (element ⇒ mass weight fraction): `Ba` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Barium()
    mat = Material("barium")
    mat.set_density(3.5)
    mat.set_mean_excitation_energy(491) # eV
    mat.add_element("Ba",1.000000)
    return mat
end

"""
    Lanthanum()

Build the predefined NIST material **LANTHANUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.154` g/cm³
- Mean excitation energy ``I``: `501` eV
- Composition (element ⇒ mass weight fraction): `La` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lanthanum()
    mat = Material("lanthanum")
    mat.set_density(6.154)
    mat.set_mean_excitation_energy(501) # eV
    mat.add_element("La",1.000000)
    return mat
end

"""
    Cerium()

Build the predefined NIST material **CERIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.657` g/cm³
- Mean excitation energy ``I``: `523` eV
- Composition (element ⇒ mass weight fraction): `Ce` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cerium()
    mat = Material("cerium")
    mat.set_density(6.657)
    mat.set_mean_excitation_energy(523) # eV
    mat.add_element("Ce",1.000000)
    return mat
end

"""
    Praseodymium()

Build the predefined NIST material **PRASEODYMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.71` g/cm³
- Mean excitation energy ``I``: `535` eV
- Composition (element ⇒ mass weight fraction): `Pr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Praseodymium()
    mat = Material("praseodymium")
    mat.set_density(6.71)
    mat.set_mean_excitation_energy(535) # eV
    mat.add_element("Pr",1.000000)
    return mat
end

"""
    Neodymium()

Build the predefined NIST material **NEODYMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.9` g/cm³
- Mean excitation energy ``I``: `546` eV
- Composition (element ⇒ mass weight fraction): `Nd` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Neodymium()
    mat = Material("neodymium")
    mat.set_density(6.9)
    mat.set_mean_excitation_energy(546) # eV
    mat.add_element("Nd",1.000000)
    return mat
end

"""
    Promethium()

Build the predefined NIST material **PROMETHIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.22` g/cm³
- Mean excitation energy ``I``: `560` eV
- Composition (element ⇒ mass weight fraction): `Pm` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Promethium()
    mat = Material("promethium")
    mat.set_density(7.22)
    mat.set_mean_excitation_energy(560) # eV
    mat.add_element("Pm",1.000000)
    return mat
end

"""
    Samarium()

Build the predefined NIST material **SAMARIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.46` g/cm³
- Mean excitation energy ``I``: `574` eV
- Composition (element ⇒ mass weight fraction): `Sm` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Samarium()
    mat = Material("samarium")
    mat.set_density(7.46)
    mat.set_mean_excitation_energy(574) # eV
    mat.add_element("Sm",1.000000)
    return mat
end

"""
    Europium()

Build the predefined NIST material **EUROPIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.243` g/cm³
- Mean excitation energy ``I``: `580` eV
- Composition (element ⇒ mass weight fraction): `Eu` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Europium()
    mat = Material("europium")
    mat.set_density(5.243)
    mat.set_mean_excitation_energy(580) # eV
    mat.add_element("Eu",1.000000)
    return mat
end

"""
    Gadolinium()

Build the predefined NIST material **GADOLINIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.9004` g/cm³
- Mean excitation energy ``I``: `591` eV
- Composition (element ⇒ mass weight fraction): `Gd` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gadolinium()
    mat = Material("gadolinium")
    mat.set_density(7.9004)
    mat.set_mean_excitation_energy(591) # eV
    mat.add_element("Gd",1.000000)
    return mat
end

"""
    Terbium()

Build the predefined NIST material **TERBIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.229` g/cm³
- Mean excitation energy ``I``: `614` eV
- Composition (element ⇒ mass weight fraction): `Tb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Terbium()
    mat = Material("terbium")
    mat.set_density(8.229)
    mat.set_mean_excitation_energy(614) # eV
    mat.add_element("Tb",1.000000)
    return mat
end

"""
    Dysprosium()

Build the predefined NIST material **DYSPROSIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.55` g/cm³
- Mean excitation energy ``I``: `628` eV
- Composition (element ⇒ mass weight fraction): `Dy` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Dysprosium()
    mat = Material("dysprosium")
    mat.set_density(8.55)
    mat.set_mean_excitation_energy(628) # eV
    mat.add_element("Dy",1.000000)
    return mat
end

"""
    Holmium()

Build the predefined NIST material **HOLMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `8.795` g/cm³
- Mean excitation energy ``I``: `650` eV
- Composition (element ⇒ mass weight fraction): `Ho` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Holmium()
    mat = Material("holmium")
    mat.set_density(8.795)
    mat.set_mean_excitation_energy(650) # eV
    mat.add_element("Ho",1.000000)
    return mat
end

"""
    Erbium()

Build the predefined NIST material **ERBIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.066` g/cm³
- Mean excitation energy ``I``: `658` eV
- Composition (element ⇒ mass weight fraction): `Er` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Erbium()
    mat = Material("erbium")
    mat.set_density(9.066)
    mat.set_mean_excitation_energy(658) # eV
    mat.add_element("Er",1.000000)
    return mat
end

"""
    Thulium()

Build the predefined NIST material **THULIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.321` g/cm³
- Mean excitation energy ``I``: `674` eV
- Composition (element ⇒ mass weight fraction): `Tm` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Thulium()
    mat = Material("thulium")
    mat.set_density(9.321)
    mat.set_mean_excitation_energy(674) # eV
    mat.add_element("Tm",1.000000)
    return mat
end

"""
    Ytterbium()

Build the predefined NIST material **YTTERBIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.73` g/cm³
- Mean excitation energy ``I``: `684` eV
- Composition (element ⇒ mass weight fraction): `Yb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ytterbium()
    mat = Material("ytterbium")
    mat.set_density(6.73)
    mat.set_mean_excitation_energy(684) # eV
    mat.add_element("Yb",1.000000)
    return mat
end

"""
    Lutetium()

Build the predefined NIST material **LUTETIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.84` g/cm³
- Mean excitation energy ``I``: `694` eV
- Composition (element ⇒ mass weight fraction): `Lu` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lutetium()
    mat = Material("lutetium")
    mat.set_density(9.84)
    mat.set_mean_excitation_energy(694) # eV
    mat.add_element("Lu",1.000000)
    return mat
end

"""
    Hafnium()

Build the predefined NIST material **HAFNIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `13.31` g/cm³
- Mean excitation energy ``I``: `705` eV
- Composition (element ⇒ mass weight fraction): `Hf` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Hafnium()
    mat = Material("hafnium")
    mat.set_density(13.31)
    mat.set_mean_excitation_energy(705) # eV
    mat.add_element("Hf",1.000000)
    return mat
end

"""
    Tantalum()

Build the predefined NIST material **TANTALUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `16.654` g/cm³
- Mean excitation energy ``I``: `718` eV
- Composition (element ⇒ mass weight fraction): `Ta` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tantalum()
    mat = Material("tantalum")
    mat.set_density(16.654)
    mat.set_mean_excitation_energy(718) # eV
    mat.add_element("Ta",1.000000)
    return mat
end

"""
    Tungsten()

Build the predefined NIST material **TUNGSTEN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `19.3` g/cm³
- Mean excitation energy ``I``: `727` eV
- Composition (element ⇒ mass weight fraction): `W` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tungsten()
    mat = Material("tungsten")
    mat.set_density(19.3)
    mat.set_mean_excitation_energy(727) # eV
    mat.add_element("W",1.000000)
    return mat
end

"""
    Rhenium()

Build the predefined NIST material **RHENIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `21.02` g/cm³
- Mean excitation energy ``I``: `736` eV
- Composition (element ⇒ mass weight fraction): `Re` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rhenium()
    mat = Material("rhenium")
    mat.set_density(21.02)
    mat.set_mean_excitation_energy(736) # eV
    mat.add_element("Re",1.000000)
    return mat
end

"""
    Osmium()

Build the predefined NIST material **OSMIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `22.57` g/cm³
- Mean excitation energy ``I``: `746` eV
- Composition (element ⇒ mass weight fraction): `Os` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Osmium()
    mat = Material("osmium")
    mat.set_density(22.57)
    mat.set_mean_excitation_energy(746) # eV
    mat.add_element("Os",1.000000)
    return mat
end

"""
    Iridium()

Build the predefined NIST material **IRIDIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `22.42` g/cm³
- Mean excitation energy ``I``: `757` eV
- Composition (element ⇒ mass weight fraction): `Ir` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Iridium()
    mat = Material("iridium")
    mat.set_density(22.42)
    mat.set_mean_excitation_energy(757) # eV
    mat.add_element("Ir",1.000000)
    return mat
end

"""
    Platinum()

Build the predefined NIST material **PLATINUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `21.45` g/cm³
- Mean excitation energy ``I``: `790` eV
- Composition (element ⇒ mass weight fraction): `Pt` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Platinum()
    mat = Material("platinum")
    mat.set_density(21.45)
    mat.set_mean_excitation_energy(790) # eV
    mat.add_element("Pt",1.000000)
    return mat
end

"""
    Gold()

Build the predefined NIST material **GOLD** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `19.32` g/cm³
- Mean excitation energy ``I``: `790` eV
- Composition (element ⇒ mass weight fraction): `Au` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gold()
    mat = Material("gold")
    mat.set_density(19.32)
    mat.set_mean_excitation_energy(790) # eV
    mat.add_element("Au",1.000000)
    return mat
end

"""
    Mercury()

Build the predefined NIST material **MERCURY** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `13.546` g/cm³
- Mean excitation energy ``I``: `800` eV
- Composition (element ⇒ mass weight fraction): `Hg` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Mercury()
    mat = Material("mercury")
    mat.set_density(13.546)
    mat.set_mean_excitation_energy(800) # eV
    mat.add_element("Hg",1.000000)
    return mat
end

"""
    Thallium()

Build the predefined NIST material **THALLIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.72` g/cm³
- Mean excitation energy ``I``: `810` eV
- Composition (element ⇒ mass weight fraction): `Tl` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Thallium()
    mat = Material("thallium")
    mat.set_density(11.72)
    mat.set_mean_excitation_energy(810) # eV
    mat.add_element("Tl",1.000000)
    return mat
end

"""
    Lead()

Build the predefined NIST material **LEAD** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.35` g/cm³
- Mean excitation energy ``I``: `823` eV
- Composition (element ⇒ mass weight fraction): `Pb` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lead()
    mat = Material("lead")
    mat.set_density(11.35)
    mat.set_mean_excitation_energy(823) # eV
    mat.add_element("Pb",1.000000)
    return mat
end

"""
    Bismuth()

Build the predefined NIST material **BISMUTH** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.747` g/cm³
- Mean excitation energy ``I``: `823` eV
- Composition (element ⇒ mass weight fraction): `Bi` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bismuth()
    mat = Material("bismuth")
    mat.set_density(9.747)
    mat.set_mean_excitation_energy(823) # eV
    mat.add_element("Bi",1.000000)
    return mat
end

"""
    Polonium()

Build the predefined NIST material **POLONIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.32` g/cm³
- Mean excitation energy ``I``: `830` eV
- Composition (element ⇒ mass weight fraction): `Po` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polonium()
    mat = Material("polonium")
    mat.set_density(9.32)
    mat.set_mean_excitation_energy(830) # eV
    mat.add_element("Po",1.000000)
    return mat
end

"""
    Astatine()

Build the predefined NIST material **ASTATINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.32` g/cm³
- Mean excitation energy ``I``: `825` eV
- Composition (element ⇒ mass weight fraction): `At` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Astatine()
    mat = Material("astatine")
    mat.set_density(9.32)
    mat.set_mean_excitation_energy(825) # eV
    mat.add_element("At",1.000000)
    return mat
end

"""
    Radon()

Build the predefined NIST material **RADON** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00906618` g/cm³
- Mean excitation energy ``I``: `794` eV
- Composition (element ⇒ mass weight fraction): `Rn` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Radon()
    mat = Material("radon")
    mat.set_density(0.00906618)
    mat.set_mean_excitation_energy(794) # eV
    mat.add_element("Rn",1.000000)
    return mat
end

"""
    Francium()

Build the predefined NIST material **FRANCIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1` g/cm³
- Mean excitation energy ``I``: `827` eV
- Composition (element ⇒ mass weight fraction): `Fr` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Francium()
    mat = Material("francium")
    mat.set_density(1)
    mat.set_mean_excitation_energy(827) # eV
    mat.add_element("Fr",1.000000)
    return mat
end

"""
    Radium()

Build the predefined NIST material **RADIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5` g/cm³
- Mean excitation energy ``I``: `826` eV
- Composition (element ⇒ mass weight fraction): `Ra` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Radium()
    mat = Material("radium")
    mat.set_density(5)
    mat.set_mean_excitation_energy(826) # eV
    mat.add_element("Ra",1.000000)
    return mat
end

"""
    Actinium()

Build the predefined NIST material **ACTINIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `10.07` g/cm³
- Mean excitation energy ``I``: `841` eV
- Composition (element ⇒ mass weight fraction): `Ac` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Actinium()
    mat = Material("actinium")
    mat.set_density(10.07)
    mat.set_mean_excitation_energy(841) # eV
    mat.add_element("Ac",1.000000)
    return mat
end

"""
    Thorium()

Build the predefined NIST material **THORIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.72` g/cm³
- Mean excitation energy ``I``: `847` eV
- Composition (element ⇒ mass weight fraction): `Th` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Thorium()
    mat = Material("thorium")
    mat.set_density(11.72)
    mat.set_mean_excitation_energy(847) # eV
    mat.add_element("Th",1.000000)
    return mat
end

"""
    Protactinium()

Build the predefined NIST material **PROTACTINIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `15.37` g/cm³
- Mean excitation energy ``I``: `878` eV
- Composition (element ⇒ mass weight fraction): `Pa` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Protactinium()
    mat = Material("protactinium")
    mat.set_density(15.37)
    mat.set_mean_excitation_energy(878) # eV
    mat.add_element("Pa",1.000000)
    return mat
end

"""
    Uranium()

Build the predefined NIST material **URANIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `18.95` g/cm³
- Mean excitation energy ``I``: `890` eV
- Composition (element ⇒ mass weight fraction): `U` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Uranium()
    mat = Material("uranium")
    mat.set_density(18.95)
    mat.set_mean_excitation_energy(890) # eV
    mat.add_element("U",1.000000)
    return mat
end

"""
    Neptunium()

Build the predefined NIST material **NEPTUNIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `20.25` g/cm³
- Mean excitation energy ``I``: `902` eV
- Composition (element ⇒ mass weight fraction): `Np` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Neptunium()
    mat = Material("neptunium")
    mat.set_density(20.25)
    mat.set_mean_excitation_energy(902) # eV
    mat.add_element("Np",1.000000)
    return mat
end

"""
    Plutonium()

Build the predefined NIST material **PLUTONIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `19.84` g/cm³
- Mean excitation energy ``I``: `921` eV
- Composition (element ⇒ mass weight fraction): `Pu` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Plutonium()
    mat = Material("plutonium")
    mat.set_density(19.84)
    mat.set_mean_excitation_energy(921) # eV
    mat.add_element("Pu",1.000000)
    return mat
end

"""
    Americium()

Build the predefined NIST material **AMERICIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `13.67` g/cm³
- Mean excitation energy ``I``: `934` eV
- Composition (element ⇒ mass weight fraction): `Am` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Americium()
    mat = Material("americium")
    mat.set_density(13.67)
    mat.set_mean_excitation_energy(934) # eV
    mat.add_element("Am",1.000000)
    return mat
end

"""
    Curium()

Build the predefined NIST material **CURIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `13.51` g/cm³
- Mean excitation energy ``I``: `939` eV
- Composition (element ⇒ mass weight fraction): `Cm` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Curium()
    mat = Material("curium")
    mat.set_density(13.51)
    mat.set_mean_excitation_energy(939) # eV
    mat.add_element("Cm",1.000000)
    return mat
end

"""
    Berkelium()

Build the predefined NIST material **BERKELIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `14` g/cm³
- Mean excitation energy ``I``: `952` eV
- Composition (element ⇒ mass weight fraction): `Bk` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Berkelium()
    mat = Material("berkelium")
    mat.set_density(14)
    mat.set_mean_excitation_energy(952) # eV
    mat.add_element("Bk",1.000000)
    return mat
end

"""
    Californium()

Build the predefined NIST material **CALIFORNIUM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `10` g/cm³
- Mean excitation energy ``I``: `966` eV
- Composition (element ⇒ mass weight fraction): `Cf` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Californium()
    mat = Material("californium")
    mat.set_density(10)
    mat.set_mean_excitation_energy(966) # eV
    mat.add_element("Cf",1.000000)
    return mat
end

"""
    A_150_Tissue_Equivalent_Plastic()

Build the predefined NIST material **A-150 TISSUE-EQUIVALENT PLASTIC** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.127` g/cm³
- Mean excitation energy ``I``: `65.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101327`
    - `C` ⇒ `0.775501`
    - `N` ⇒ `0.035057`
    - `O` ⇒ `0.052316`
    - `F` ⇒ `0.017422`
    - `Ca` ⇒ `0.018378`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function A_150_Tissue_Equivalent_Plastic()
    mat = Material("a_150_tissue_equivalent_plastic")
    mat.set_density(1.127)
    mat.set_mean_excitation_energy(65.1) # eV
    mat.add_element("H",0.101327)
    mat.add_element("C",0.775501)
    mat.add_element("N",0.035057)
    mat.add_element("O",0.052316)
    mat.add_element("F",0.017422)
    mat.add_element("Ca",0.018378)
    return mat
end

"""
    Acetone()

Build the predefined NIST material **ACETONE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.7899` g/cm³
- Mean excitation energy ``I``: `64.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.104122`
    - `C` ⇒ `0.620405`
    - `O` ⇒ `0.275473`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Acetone()
    mat = Material("acetone")
    mat.set_density(0.7899)
    mat.set_mean_excitation_energy(64.2) # eV
    mat.add_element("H",0.104122)
    mat.add_element("C",0.620405)
    mat.add_element("O",0.275473)
    return mat
end

"""
    Acetylene()

Build the predefined NIST material **ACETYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.0010967` g/cm³
- Mean excitation energy ``I``: `58.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.077418`
    - `C` ⇒ `0.922582`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Acetylene()
    mat = Material("acetylene")
    mat.set_density(0.0010967)
    mat.set_mean_excitation_energy(58.2) # eV
    mat.add_element("H",0.077418)
    mat.add_element("C",0.922582)
    return mat
end

"""
    Adenine()

Build the predefined NIST material **ADENINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.35` g/cm³
- Mean excitation energy ``I``: `71.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.037294`
    - `C` ⇒ `0.444430`
    - `N` ⇒ `0.518275`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Adenine()
    mat = Material("adenine")
    mat.set_density(1.35)
    mat.set_mean_excitation_energy(71.4) # eV
    mat.add_element("H",0.037294)
    mat.add_element("C",0.444430)
    mat.add_element("N",0.518275)
    return mat
end

"""
    Adipose_Tissue_Icrp()

Build the predefined NIST material **ADIPOSE TISSUE (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.92` g/cm³
- Mean excitation energy ``I``: `63.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.119477`
    - `C` ⇒ `0.637240`
    - `N` ⇒ `0.007970`
    - `O` ⇒ `0.232333`
    - `Na` ⇒ `0.000500`
    - `Mg` ⇒ `0.000020`
    - `P` ⇒ `0.000160`
    - `S` ⇒ `0.000730`
    - `Cl` ⇒ `0.001190`
    - `K` ⇒ `0.000320`
    - `Ca` ⇒ `0.000020`
    - `Fe` ⇒ `0.000020`
    - `Zn` ⇒ `0.000020`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Adipose_Tissue_Icrp()
    mat = Material("adipose_tissue_icrp")
    mat.set_density(0.92)
    mat.set_mean_excitation_energy(63.2) # eV
    mat.add_element("H",0.119477)
    mat.add_element("C",0.637240)
    mat.add_element("N",0.007970)
    mat.add_element("O",0.232333)
    mat.add_element("Na",0.000500)
    mat.add_element("Mg",0.000020)
    mat.add_element("P",0.000160)
    mat.add_element("S",0.000730)
    mat.add_element("Cl",0.001190)
    mat.add_element("K",0.000320)
    mat.add_element("Ca",0.000020)
    mat.add_element("Fe",0.000020)
    mat.add_element("Zn",0.000020)
    return mat
end

"""
    Air_Dry_Near_Sea_Level()

Build the predefined NIST material **AIR, DRY (NEAR SEA LEVEL)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00120479` g/cm³
- Mean excitation energy ``I``: `85.7` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.000124`
    - `N` ⇒ `0.755267`
    - `O` ⇒ `0.231781`
    - `Ar` ⇒ `0.012827`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Air_Dry_Near_Sea_Level()
    mat = Material("air_dry_near_sea_level")
    mat.set_density(0.00120479)
    mat.set_mean_excitation_energy(85.7) # eV
    mat.add_element("C",0.000124)
    mat.add_element("N",0.755267)
    mat.add_element("O",0.231781)
    mat.add_element("Ar",0.012827)
    return mat
end

"""
    Alanine()

Build the predefined NIST material **ALANINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.42` g/cm³
- Mean excitation energy ``I``: `71.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.079190`
    - `C` ⇒ `0.404439`
    - `N` ⇒ `0.157213`
    - `O` ⇒ `0.359159`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Alanine()
    mat = Material("alanine")
    mat.set_density(1.42)
    mat.set_mean_excitation_energy(71.9) # eV
    mat.add_element("H",0.079190)
    mat.add_element("C",0.404439)
    mat.add_element("N",0.157213)
    mat.add_element("O",0.359159)
    return mat
end

"""
    Aluminum_Oxide()

Build the predefined NIST material **ALUMINUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.97` g/cm³
- Mean excitation energy ``I``: `145.2` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.470749`
    - `Al` ⇒ `0.529251`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Aluminum_Oxide()
    mat = Material("aluminum_oxide")
    mat.set_density(3.97)
    mat.set_mean_excitation_energy(145.2) # eV
    mat.add_element("O",0.470749)
    mat.add_element("Al",0.529251)
    return mat
end

"""
    Amber()

Build the predefined NIST material **AMBER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.1` g/cm³
- Mean excitation energy ``I``: `63.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.105930`
    - `C` ⇒ `0.788973`
    - `O` ⇒ `0.105096`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Amber()
    mat = Material("amber")
    mat.set_density(1.1)
    mat.set_mean_excitation_energy(63.2) # eV
    mat.add_element("H",0.105930)
    mat.add_element("C",0.788973)
    mat.add_element("O",0.105096)
    return mat
end

"""
    Ammonia()

Build the predefined NIST material **AMMONIA** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.000826019` g/cm³
- Mean excitation energy ``I``: `53.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.177547`
    - `N` ⇒ `0.822453`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ammonia()
    mat = Material("ammonia")
    mat.set_density(0.000826019)
    mat.set_mean_excitation_energy(53.7) # eV
    mat.add_element("H",0.177547)
    mat.add_element("N",0.822453)
    return mat
end

"""
    Aniline()

Build the predefined NIST material **ANILINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.0235` g/cm³
- Mean excitation energy ``I``: `66.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.075759`
    - `C` ⇒ `0.773838`
    - `N` ⇒ `0.150403`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Aniline()
    mat = Material("aniline")
    mat.set_density(1.0235)
    mat.set_mean_excitation_energy(66.2) # eV
    mat.add_element("H",0.075759)
    mat.add_element("C",0.773838)
    mat.add_element("N",0.150403)
    return mat
end

"""
    Anthracene()

Build the predefined NIST material **ANTHRACENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.283` g/cm³
- Mean excitation energy ``I``: `69.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.056550`
    - `C` ⇒ `0.943450`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Anthracene()
    mat = Material("anthracene")
    mat.set_density(1.283)
    mat.set_mean_excitation_energy(69.5) # eV
    mat.add_element("H",0.056550)
    mat.add_element("C",0.943450)
    return mat
end

"""
    B_100_Bone_Equivalent_Plastic()

Build the predefined NIST material **B-100 BONE-EQUIVALENT PLASTIC** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.45` g/cm³
- Mean excitation energy ``I``: `85.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.065471`
    - `C` ⇒ `0.536945`
    - `N` ⇒ `0.021500`
    - `O` ⇒ `0.032085`
    - `F` ⇒ `0.167411`
    - `Ca` ⇒ `0.176589`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function B_100_Bone_Equivalent_Plastic()
    mat = Material("b_100_bone_equivalent_plastic")
    mat.set_density(1.45)
    mat.set_mean_excitation_energy(85.9) # eV
    mat.add_element("H",0.065471)
    mat.add_element("C",0.536945)
    mat.add_element("N",0.021500)
    mat.add_element("O",0.032085)
    mat.add_element("F",0.167411)
    mat.add_element("Ca",0.176589)
    return mat
end

"""
    Bakelite()

Build the predefined NIST material **BAKELITE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.25` g/cm³
- Mean excitation energy ``I``: `72.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.057441`
    - `C` ⇒ `0.774591`
    - `O` ⇒ `0.167968`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bakelite()
    mat = Material("bakelite")
    mat.set_density(1.25)
    mat.set_mean_excitation_energy(72.4) # eV
    mat.add_element("H",0.057441)
    mat.add_element("C",0.774591)
    mat.add_element("O",0.167968)
    return mat
end

"""
    Barium_Fluoride()

Build the predefined NIST material **BARIUM FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.89` g/cm³
- Mean excitation energy ``I``: `375.9` eV
- Composition (element ⇒ mass weight fraction):
    - `F` ⇒ `0.216720`
    - `Ba` ⇒ `0.783280`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Barium_Fluoride()
    mat = Material("barium_fluoride")
    mat.set_density(4.89)
    mat.set_mean_excitation_energy(375.9) # eV
    mat.add_element("F",0.216720)
    mat.add_element("Ba",0.783280)
    return mat
end

"""
    Barium_Sulfate()

Build the predefined NIST material **BARIUM SULFATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.5` g/cm³
- Mean excitation energy ``I``: `285.7` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.274212`
    - `S` ⇒ `0.137368`
    - `Ba` ⇒ `0.588420`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Barium_Sulfate()
    mat = Material("barium_sulfate")
    mat.set_density(4.5)
    mat.set_mean_excitation_energy(285.7) # eV
    mat.add_element("O",0.274212)
    mat.add_element("S",0.137368)
    mat.add_element("Ba",0.588420)
    return mat
end

"""
    Benzene()

Build the predefined NIST material **BENZENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.87865` g/cm³
- Mean excitation energy ``I``: `63.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.077418`
    - `C` ⇒ `0.922582`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Benzene()
    mat = Material("benzene")
    mat.set_density(0.87865)
    mat.set_mean_excitation_energy(63.4) # eV
    mat.add_element("H",0.077418)
    mat.add_element("C",0.922582)
    return mat
end

"""
    Beryllium_Oxide()

Build the predefined NIST material **BERYLLIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.01` g/cm³
- Mean excitation energy ``I``: `93.2` eV
- Composition (element ⇒ mass weight fraction):
    - `Be` ⇒ `0.360320`
    - `O` ⇒ `0.639680`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Beryllium_Oxide()
    mat = Material("beryllium_oxide")
    mat.set_density(3.01)
    mat.set_mean_excitation_energy(93.2) # eV
    mat.add_element("Be",0.360320)
    mat.add_element("O",0.639680)
    return mat
end

"""
    Bismuth_Germanium_Oxide()

Build the predefined NIST material **BISMUTH GERMANIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.13` g/cm³
- Mean excitation energy ``I``: `534.1` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.154126`
    - `Ge` ⇒ `0.174820`
    - `Bi` ⇒ `0.671054`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bismuth_Germanium_Oxide()
    mat = Material("bismuth_germanium_oxide")
    mat.set_density(7.13)
    mat.set_mean_excitation_energy(534.1) # eV
    mat.add_element("O",0.154126)
    mat.add_element("Ge",0.174820)
    mat.add_element("Bi",0.671054)
    return mat
end

"""
    Blood_Icrp()

Build the predefined NIST material **BLOOD (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.06` g/cm³
- Mean excitation energy ``I``: `75.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101866`
    - `C` ⇒ `0.100020`
    - `N` ⇒ `0.029640`
    - `O` ⇒ `0.759414`
    - `Na` ⇒ `0.001850`
    - `Mg` ⇒ `0.000040`
    - `Si` ⇒ `0.000030`
    - `P` ⇒ `0.000350`
    - `S` ⇒ `0.001850`
    - `Cl` ⇒ `0.002780`
    - `K` ⇒ `0.001630`
    - `Ca` ⇒ `0.000060`
    - `Fe` ⇒ `0.000460`
    - `Zn` ⇒ `0.000010`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Blood_Icrp()
    mat = Material("blood_icrp")
    mat.set_density(1.06)
    mat.set_mean_excitation_energy(75.2) # eV
    mat.add_element("H",0.101866)
    mat.add_element("C",0.100020)
    mat.add_element("N",0.029640)
    mat.add_element("O",0.759414)
    mat.add_element("Na",0.001850)
    mat.add_element("Mg",0.000040)
    mat.add_element("Si",0.000030)
    mat.add_element("P",0.000350)
    mat.add_element("S",0.001850)
    mat.add_element("Cl",0.002780)
    mat.add_element("K",0.001630)
    mat.add_element("Ca",0.000060)
    mat.add_element("Fe",0.000460)
    mat.add_element("Zn",0.000010)
    return mat
end

"""
    Bone_Compact_Icru()

Build the predefined NIST material **BONE, COMPACT (ICRU)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.85` g/cm³
- Mean excitation energy ``I``: `91.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.063984`
    - `C` ⇒ `0.278000`
    - `N` ⇒ `0.027000`
    - `O` ⇒ `0.410016`
    - `Mg` ⇒ `0.002000`
    - `P` ⇒ `0.070000`
    - `S` ⇒ `0.002000`
    - `Ca` ⇒ `0.147000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bone_Compact_Icru()
    mat = Material("bone_compact_icru")
    mat.set_density(1.85)
    mat.set_mean_excitation_energy(91.9) # eV
    mat.add_element("H",0.063984)
    mat.add_element("C",0.278000)
    mat.add_element("N",0.027000)
    mat.add_element("O",0.410016)
    mat.add_element("Mg",0.002000)
    mat.add_element("P",0.070000)
    mat.add_element("S",0.002000)
    mat.add_element("Ca",0.147000)
    return mat
end

"""
    Bone_Cortical_Icrp()

Build the predefined NIST material **BONE, CORTICAL (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.85` g/cm³
- Mean excitation energy ``I``: `106.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.047234`
    - `C` ⇒ `0.144330`
    - `N` ⇒ `0.041990`
    - `O` ⇒ `0.446096`
    - `Mg` ⇒ `0.002200`
    - `P` ⇒ `0.104970`
    - `S` ⇒ `0.003150`
    - `Ca` ⇒ `0.209930`
    - `Zn` ⇒ `0.000100`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Bone_Cortical_Icrp()
    mat = Material("bone_cortical_icrp")
    mat.set_density(1.85)
    mat.set_mean_excitation_energy(106.4) # eV
    mat.add_element("H",0.047234)
    mat.add_element("C",0.144330)
    mat.add_element("N",0.041990)
    mat.add_element("O",0.446096)
    mat.add_element("Mg",0.002200)
    mat.add_element("P",0.104970)
    mat.add_element("S",0.003150)
    mat.add_element("Ca",0.209930)
    mat.add_element("Zn",0.000100)
    return mat
end

"""
    Boron_Carbide()

Build the predefined NIST material **BORON CARBIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.52` g/cm³
- Mean excitation energy ``I``: `84.7` eV
- Composition (element ⇒ mass weight fraction):
    - `B` ⇒ `0.782610`
    - `C` ⇒ `0.217390`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Boron_Carbide()
    mat = Material("boron_carbide")
    mat.set_density(2.52)
    mat.set_mean_excitation_energy(84.7) # eV
    mat.add_element("B",0.782610)
    mat.add_element("C",0.217390)
    return mat
end

"""
    Boron_Oxide()

Build the predefined NIST material **BORON OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.812` g/cm³
- Mean excitation energy ``I``: `99.6` eV
- Composition (element ⇒ mass weight fraction):
    - `B` ⇒ `0.310551`
    - `O` ⇒ `0.689449`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Boron_Oxide()
    mat = Material("boron_oxide")
    mat.set_density(1.812)
    mat.set_mean_excitation_energy(99.6) # eV
    mat.add_element("B",0.310551)
    mat.add_element("O",0.689449)
    return mat
end

"""
    Brain_Icrp()

Build the predefined NIST material **BRAIN (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.03` g/cm³
- Mean excitation energy ``I``: `73.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.110667`
    - `C` ⇒ `0.125420`
    - `N` ⇒ `0.013280`
    - `O` ⇒ `0.737723`
    - `Na` ⇒ `0.001840`
    - `Mg` ⇒ `0.000150`
    - `P` ⇒ `0.003540`
    - `S` ⇒ `0.001770`
    - `Cl` ⇒ `0.002360`
    - `K` ⇒ `0.003100`
    - `Ca` ⇒ `0.000090`
    - `Fe` ⇒ `0.000050`
    - `Zn` ⇒ `0.000010`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Brain_Icrp()
    mat = Material("brain_icrp")
    mat.set_density(1.03)
    mat.set_mean_excitation_energy(73.3) # eV
    mat.add_element("H",0.110667)
    mat.add_element("C",0.125420)
    mat.add_element("N",0.013280)
    mat.add_element("O",0.737723)
    mat.add_element("Na",0.001840)
    mat.add_element("Mg",0.000150)
    mat.add_element("P",0.003540)
    mat.add_element("S",0.001770)
    mat.add_element("Cl",0.002360)
    mat.add_element("K",0.003100)
    mat.add_element("Ca",0.000090)
    mat.add_element("Fe",0.000050)
    mat.add_element("Zn",0.000010)
    return mat
end

"""
    Butane()

Build the predefined NIST material **BUTANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00249343` g/cm³
- Mean excitation energy ``I``: `48.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.173408`
    - `C` ⇒ `0.826592`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Butane()
    mat = Material("butane")
    mat.set_density(0.00249343)
    mat.set_mean_excitation_energy(48.3) # eV
    mat.add_element("H",0.173408)
    mat.add_element("C",0.826592)
    return mat
end

"""
    N_Butyl_Alcohol()

Build the predefined NIST material **N-BUTYL ALCOHOL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.8098` g/cm³
- Mean excitation energy ``I``: `59.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.135978`
    - `C` ⇒ `0.648171`
    - `O` ⇒ `0.215851`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_Butyl_Alcohol()
    mat = Material("n_butyl_alcohol")
    mat.set_density(0.8098)
    mat.set_mean_excitation_energy(59.9) # eV
    mat.add_element("H",0.135978)
    mat.add_element("C",0.648171)
    mat.add_element("O",0.215851)
    return mat
end

"""
    C_552_Air_Equivalent_Plastic()

Build the predefined NIST material **C-552 AIR-EQUIVALENT PLASTIC** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.76` g/cm³
- Mean excitation energy ``I``: `86.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.024680`
    - `C` ⇒ `0.501610`
    - `O` ⇒ `0.004527`
    - `F` ⇒ `0.465209`
    - `Si` ⇒ `0.003973`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function C_552_Air_Equivalent_Plastic()
    mat = Material("c_552_air_equivalent_plastic")
    mat.set_density(1.76)
    mat.set_mean_excitation_energy(86.8) # eV
    mat.add_element("H",0.024680)
    mat.add_element("C",0.501610)
    mat.add_element("O",0.004527)
    mat.add_element("F",0.465209)
    mat.add_element("Si",0.003973)
    return mat
end

"""
    Cadmium_Telluride()

Build the predefined NIST material **CADMIUM TELLURIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.2` g/cm³
- Mean excitation energy ``I``: `539.3` eV
- Composition (element ⇒ mass weight fraction):
    - `Cd` ⇒ `0.468355`
    - `Te` ⇒ `0.531645`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cadmium_Telluride()
    mat = Material("cadmium_telluride")
    mat.set_density(6.2)
    mat.set_mean_excitation_energy(539.3) # eV
    mat.add_element("Cd",0.468355)
    mat.add_element("Te",0.531645)
    return mat
end

"""
    Cadmium_Tungstate()

Build the predefined NIST material **CADMIUM TUNGSTATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.9` g/cm³
- Mean excitation energy ``I``: `468.3` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.177644`
    - `Cd` ⇒ `0.312027`
    - `W` ⇒ `0.510329`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cadmium_Tungstate()
    mat = Material("cadmium_tungstate")
    mat.set_density(7.9)
    mat.set_mean_excitation_energy(468.3) # eV
    mat.add_element("O",0.177644)
    mat.add_element("Cd",0.312027)
    mat.add_element("W",0.510329)
    return mat
end

"""
    Calcium_Carbonate()

Build the predefined NIST material **CALCIUM CARBONATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.8` g/cm³
- Mean excitation energy ``I``: `136.4` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.120003`
    - `O` ⇒ `0.479554`
    - `Ca` ⇒ `0.400443`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium_Carbonate()
    mat = Material("calcium_carbonate")
    mat.set_density(2.8)
    mat.set_mean_excitation_energy(136.4) # eV
    mat.add_element("C",0.120003)
    mat.add_element("O",0.479554)
    mat.add_element("Ca",0.400443)
    return mat
end

"""
    Calcium_Fluoride()

Build the predefined NIST material **CALCIUM FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.18` g/cm³
- Mean excitation energy ``I``: `166` eV
- Composition (element ⇒ mass weight fraction):
    - `F` ⇒ `0.486659`
    - `Ca` ⇒ `0.513341`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium_Fluoride()
    mat = Material("calcium_fluoride")
    mat.set_density(3.18)
    mat.set_mean_excitation_energy(166) # eV
    mat.add_element("F",0.486659)
    mat.add_element("Ca",0.513341)
    return mat
end

"""
    Calcium_Oxide()

Build the predefined NIST material **CALCIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.3` g/cm³
- Mean excitation energy ``I``: `176.1` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.285299`
    - `Ca` ⇒ `0.714701`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium_Oxide()
    mat = Material("calcium_oxide")
    mat.set_density(3.3)
    mat.set_mean_excitation_energy(176.1) # eV
    mat.add_element("O",0.285299)
    mat.add_element("Ca",0.714701)
    return mat
end

"""
    Calcium_Sulfate()

Build the predefined NIST material **CALCIUM SULFATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.96` g/cm³
- Mean excitation energy ``I``: `152.3` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.470095`
    - `S` ⇒ `0.235497`
    - `Ca` ⇒ `0.294408`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium_Sulfate()
    mat = Material("calcium_sulfate")
    mat.set_density(2.96)
    mat.set_mean_excitation_energy(152.3) # eV
    mat.add_element("O",0.470095)
    mat.add_element("S",0.235497)
    mat.add_element("Ca",0.294408)
    return mat
end

"""
    Calcium_Tungstate()

Build the predefined NIST material **CALCIUM TUNGSTATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.062` g/cm³
- Mean excitation energy ``I``: `395` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.222270`
    - `Ca` ⇒ `0.139202`
    - `W` ⇒ `0.638529`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Calcium_Tungstate()
    mat = Material("calcium_tungstate")
    mat.set_density(6.062)
    mat.set_mean_excitation_energy(395) # eV
    mat.add_element("O",0.222270)
    mat.add_element("Ca",0.139202)
    mat.add_element("W",0.638529)
    return mat
end

"""
    Carbon_Dioxide()

Build the predefined NIST material **CARBON DIOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00184212` g/cm³
- Mean excitation energy ``I``: `85` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.272916`
    - `O` ⇒ `0.727084`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Carbon_Dioxide()
    mat = Material("carbon_dioxide")
    mat.set_density(0.00184212)
    mat.set_mean_excitation_energy(85) # eV
    mat.add_element("C",0.272916)
    mat.add_element("O",0.727084)
    return mat
end

"""
    Carbon_Tetrachloride()

Build the predefined NIST material **CARBON TETRACHLORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.594` g/cm³
- Mean excitation energy ``I``: `166.3` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.078083`
    - `Cl` ⇒ `0.921917`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Carbon_Tetrachloride()
    mat = Material("carbon_tetrachloride")
    mat.set_density(1.594)
    mat.set_mean_excitation_energy(166.3) # eV
    mat.add_element("C",0.078083)
    mat.add_element("Cl",0.921917)
    return mat
end

"""
    Cellulose_Acetate_Cellophane()

Build the predefined NIST material **CELLULOSE ACETATE, CELLOPHANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.42` g/cm³
- Mean excitation energy ``I``: `77.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.062162`
    - `C` ⇒ `0.444462`
    - `O` ⇒ `0.493376`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cellulose_Acetate_Cellophane()
    mat = Material("cellulose_acetate_cellophane")
    mat.set_density(1.42)
    mat.set_mean_excitation_energy(77.6) # eV
    mat.add_element("H",0.062162)
    mat.add_element("C",0.444462)
    mat.add_element("O",0.493376)
    return mat
end

"""
    Cellulose_Acetate_Butyrate()

Build the predefined NIST material **CELLULOSE ACETATE BUTYRATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2` g/cm³
- Mean excitation energy ``I``: `74.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.067125`
    - `C` ⇒ `0.545403`
    - `O` ⇒ `0.387472`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cellulose_Acetate_Butyrate()
    mat = Material("cellulose_acetate_butyrate")
    mat.set_density(1.2)
    mat.set_mean_excitation_energy(74.6) # eV
    mat.add_element("H",0.067125)
    mat.add_element("C",0.545403)
    mat.add_element("O",0.387472)
    return mat
end

"""
    Cellulose_Nitrate()

Build the predefined NIST material **CELLULOSE NITRATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.49` g/cm³
- Mean excitation energy ``I``: `87` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.029216`
    - `C` ⇒ `0.271296`
    - `N` ⇒ `0.121276`
    - `O` ⇒ `0.578212`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cellulose_Nitrate()
    mat = Material("cellulose_nitrate")
    mat.set_density(1.49)
    mat.set_mean_excitation_energy(87) # eV
    mat.add_element("H",0.029216)
    mat.add_element("C",0.271296)
    mat.add_element("N",0.121276)
    mat.add_element("O",0.578212)
    return mat
end

"""
    Ceric_Sulfate_Dosimeter_Solution()

Build the predefined NIST material **CERIC SULFATE DOSIMETER SOLUTION** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.03` g/cm³
- Mean excitation energy ``I``: `76.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.107596`
    - `N` ⇒ `0.000800`
    - `O` ⇒ `0.874976`
    - `S` ⇒ `0.014627`
    - `Ce` ⇒ `0.002001`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ceric_Sulfate_Dosimeter_Solution()
    mat = Material("ceric_sulfate_dosimeter_solution")
    mat.set_density(1.03)
    mat.set_mean_excitation_energy(76.7) # eV
    mat.add_element("H",0.107596)
    mat.add_element("N",0.000800)
    mat.add_element("O",0.874976)
    mat.add_element("S",0.014627)
    mat.add_element("Ce",0.002001)
    return mat
end

"""
    Cesium_Fluoride()

Build the predefined NIST material **CESIUM FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.115` g/cm³
- Mean excitation energy ``I``: `440.7` eV
- Composition (element ⇒ mass weight fraction):
    - `F` ⇒ `0.125069`
    - `Cs` ⇒ `0.874931`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cesium_Fluoride()
    mat = Material("cesium_fluoride")
    mat.set_density(4.115)
    mat.set_mean_excitation_energy(440.7) # eV
    mat.add_element("F",0.125069)
    mat.add_element("Cs",0.874931)
    return mat
end

"""
    Cesium_Iodide()

Build the predefined NIST material **CESIUM IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.51` g/cm³
- Mean excitation energy ``I``: `553.1` eV
- Composition (element ⇒ mass weight fraction):
    - `I` ⇒ `0.488451`
    - `Cs` ⇒ `0.511549`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cesium_Iodide()
    mat = Material("cesium_iodide")
    mat.set_density(4.51)
    mat.set_mean_excitation_energy(553.1) # eV
    mat.add_element("I",0.488451)
    mat.add_element("Cs",0.511549)
    return mat
end

"""
    Chlorobenzene()

Build the predefined NIST material **CHLOROBENZENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.1058` g/cm³
- Mean excitation energy ``I``: `89.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.044772`
    - `C` ⇒ `0.640254`
    - `Cl` ⇒ `0.314974`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Chlorobenzene()
    mat = Material("chlorobenzene")
    mat.set_density(1.1058)
    mat.set_mean_excitation_energy(89.1) # eV
    mat.add_element("H",0.044772)
    mat.add_element("C",0.640254)
    mat.add_element("Cl",0.314974)
    return mat
end

"""
    Chloroform()

Build the predefined NIST material **CHLOROFORM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.4832` g/cm³
- Mean excitation energy ``I``: `156` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.008443`
    - `C` ⇒ `0.100613`
    - `Cl` ⇒ `0.890944`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Chloroform()
    mat = Material("chloroform")
    mat.set_density(1.4832)
    mat.set_mean_excitation_energy(156) # eV
    mat.add_element("H",0.008443)
    mat.add_element("C",0.100613)
    mat.add_element("Cl",0.890944)
    return mat
end

"""
    Concrete_Portland()

Build the predefined NIST material **CONCRETE, PORTLAND** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.3` g/cm³
- Mean excitation energy ``I``: `135.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.010000`
    - `C` ⇒ `0.001000`
    - `O` ⇒ `0.529107`
    - `Na` ⇒ `0.016000`
    - `Mg` ⇒ `0.002000`
    - `Al` ⇒ `0.033872`
    - `Si` ⇒ `0.337021`
    - `K` ⇒ `0.013000`
    - `Ca` ⇒ `0.044000`
    - `Fe` ⇒ `0.014000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Concrete_Portland()
    mat = Material("concrete_portland")
    mat.set_density(2.3)
    mat.set_mean_excitation_energy(135.2) # eV
    mat.add_element("H",0.010000)
    mat.add_element("C",0.001000)
    mat.add_element("O",0.529107)
    mat.add_element("Na",0.016000)
    mat.add_element("Mg",0.002000)
    mat.add_element("Al",0.033872)
    mat.add_element("Si",0.337021)
    mat.add_element("K",0.013000)
    mat.add_element("Ca",0.044000)
    mat.add_element("Fe",0.014000)
    return mat
end

"""
    Cyclohexane()

Build the predefined NIST material **CYCLOHEXANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.779` g/cm³
- Mean excitation energy ``I``: `56.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.143711`
    - `C` ⇒ `0.856289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Cyclohexane()
    mat = Material("cyclohexane")
    mat.set_density(0.779)
    mat.set_mean_excitation_energy(56.4) # eV
    mat.add_element("H",0.143711)
    mat.add_element("C",0.856289)
    return mat
end

"""
    M1_2_Dichlorobenzene()

Build the predefined NIST material **1,2-DICHLOROBENZENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.3048` g/cm³
- Mean excitation energy ``I``: `106.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.027425`
    - `C` ⇒ `0.490233`
    - `Cl` ⇒ `0.482342`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function M1_2_Dichlorobenzene()
    mat = Material("1_2_dichlorobenzene")
    mat.set_density(1.3048)
    mat.set_mean_excitation_energy(106.5) # eV
    mat.add_element("H",0.027425)
    mat.add_element("C",0.490233)
    mat.add_element("Cl",0.482342)
    return mat
end

"""
    Dichlorodiethyl_Ether()

Build the predefined NIST material **DICHLORODIETHYL ETHER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2199` g/cm³
- Mean excitation energy ``I``: `103.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.056381`
    - `C` ⇒ `0.335942`
    - `O` ⇒ `0.111874`
    - `Cl` ⇒ `0.495802`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Dichlorodiethyl_Ether()
    mat = Material("dichlorodiethyl_ether")
    mat.set_density(1.2199)
    mat.set_mean_excitation_energy(103.3) # eV
    mat.add_element("H",0.056381)
    mat.add_element("C",0.335942)
    mat.add_element("O",0.111874)
    mat.add_element("Cl",0.495802)
    return mat
end

"""
    M1_2_Dichloroethane()

Build the predefined NIST material **1,2-DICHLOROETHANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2351` g/cm³
- Mean excitation energy ``I``: `111.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.040740`
    - `C` ⇒ `0.242746`
    - `Cl` ⇒ `0.716515`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function M1_2_Dichloroethane()
    mat = Material("1_2_dichloroethane")
    mat.set_density(1.2351)
    mat.set_mean_excitation_energy(111.9) # eV
    mat.add_element("H",0.040740)
    mat.add_element("C",0.242746)
    mat.add_element("Cl",0.716515)
    return mat
end

"""
    Diethyl_Ether()

Build the predefined NIST material **DIETHYL ETHER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.71378` g/cm³
- Mean excitation energy ``I``: `60` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.135978`
    - `C` ⇒ `0.648171`
    - `O` ⇒ `0.215851`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Diethyl_Ether()
    mat = Material("diethyl_ether")
    mat.set_density(0.71378)
    mat.set_mean_excitation_energy(60) # eV
    mat.add_element("H",0.135978)
    mat.add_element("C",0.648171)
    mat.add_element("O",0.215851)
    return mat
end

"""
    N_N_Dimethyl_Formamide()

Build the predefined NIST material **N,N-DIMETHYL FORMAMIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.9487` g/cm³
- Mean excitation energy ``I``: `66.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.096523`
    - `C` ⇒ `0.492965`
    - `N` ⇒ `0.191625`
    - `O` ⇒ `0.218887`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_N_Dimethyl_Formamide()
    mat = Material("n_n_dimethyl_formamide")
    mat.set_density(0.9487)
    mat.set_mean_excitation_energy(66.6) # eV
    mat.add_element("H",0.096523)
    mat.add_element("C",0.492965)
    mat.add_element("N",0.191625)
    mat.add_element("O",0.218887)
    return mat
end

"""
    Dimethyl_Sulfoxide()

Build the predefined NIST material **DIMETHYL SULFOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.1014` g/cm³
- Mean excitation energy ``I``: `98.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.077403`
    - `C` ⇒ `0.307467`
    - `O` ⇒ `0.204782`
    - `S` ⇒ `0.410348`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Dimethyl_Sulfoxide()
    mat = Material("dimethyl_sulfoxide")
    mat.set_density(1.1014)
    mat.set_mean_excitation_energy(98.6) # eV
    mat.add_element("H",0.077403)
    mat.add_element("C",0.307467)
    mat.add_element("O",0.204782)
    mat.add_element("S",0.410348)
    return mat
end

"""
    Ethane()

Build the predefined NIST material **ETHANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00125324` g/cm³
- Mean excitation energy ``I``: `45.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.201115`
    - `C` ⇒ `0.798885`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ethane()
    mat = Material("ethane")
    mat.set_density(0.00125324)
    mat.set_mean_excitation_energy(45.4) # eV
    mat.add_element("H",0.201115)
    mat.add_element("C",0.798885)
    return mat
end

"""
    Ethyl_Alcohol()

Build the predefined NIST material **ETHYL ALCOHOL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.7893` g/cm³
- Mean excitation energy ``I``: `62.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.131269`
    - `C` ⇒ `0.521438`
    - `O` ⇒ `0.347294`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ethyl_Alcohol()
    mat = Material("ethyl_alcohol")
    mat.set_density(0.7893)
    mat.set_mean_excitation_energy(62.9) # eV
    mat.add_element("H",0.131269)
    mat.add_element("C",0.521438)
    mat.add_element("O",0.347294)
    return mat
end

"""
    Ethyl_Cellulose()

Build the predefined NIST material **ETHYL CELLULOSE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.13` g/cm³
- Mean excitation energy ``I``: `69.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.090027`
    - `C` ⇒ `0.585182`
    - `O` ⇒ `0.324791`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ethyl_Cellulose()
    mat = Material("ethyl_cellulose")
    mat.set_density(1.13)
    mat.set_mean_excitation_energy(69.3) # eV
    mat.add_element("H",0.090027)
    mat.add_element("C",0.585182)
    mat.add_element("O",0.324791)
    return mat
end

"""
    Ethylene()

Build the predefined NIST material **ETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00117497` g/cm³
- Mean excitation energy ``I``: `50.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.143711`
    - `C` ⇒ `0.856289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ethylene()
    mat = Material("ethylene")
    mat.set_density(0.00117497)
    mat.set_mean_excitation_energy(50.7) # eV
    mat.add_element("H",0.143711)
    mat.add_element("C",0.856289)
    return mat
end

"""
    Eye_Lens_Icrp()

Build the predefined NIST material **EYE LENS (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.1` g/cm³
- Mean excitation energy ``I``: `73.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.099269`
    - `C` ⇒ `0.193710`
    - `N` ⇒ `0.053270`
    - `O` ⇒ `0.653751`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Eye_Lens_Icrp()
    mat = Material("eye_lens_icrp")
    mat.set_density(1.1)
    mat.set_mean_excitation_energy(73.3) # eV
    mat.add_element("H",0.099269)
    mat.add_element("C",0.193710)
    mat.add_element("N",0.053270)
    mat.add_element("O",0.653751)
    return mat
end

"""
    Ferric_Oxide()

Build the predefined NIST material **FERRIC OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.2` g/cm³
- Mean excitation energy ``I``: `227.3` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.300567`
    - `Fe` ⇒ `0.699433`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ferric_Oxide()
    mat = Material("ferric_oxide")
    mat.set_density(5.2)
    mat.set_mean_excitation_energy(227.3) # eV
    mat.add_element("O",0.300567)
    mat.add_element("Fe",0.699433)
    return mat
end

"""
    Ferroboride()

Build the predefined NIST material **FERROBORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.15` g/cm³
- Mean excitation energy ``I``: `261` eV
- Composition (element ⇒ mass weight fraction):
    - `B` ⇒ `0.162174`
    - `Fe` ⇒ `0.837826`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ferroboride()
    mat = Material("ferroboride")
    mat.set_density(7.15)
    mat.set_mean_excitation_energy(261) # eV
    mat.add_element("B",0.162174)
    mat.add_element("Fe",0.837826)
    return mat
end

"""
    Ferrous_Oxide()

Build the predefined NIST material **FERROUS OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.7` g/cm³
- Mean excitation energy ``I``: `248.6` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.222689`
    - `Fe` ⇒ `0.777311`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ferrous_Oxide()
    mat = Material("ferrous_oxide")
    mat.set_density(5.7)
    mat.set_mean_excitation_energy(248.6) # eV
    mat.add_element("O",0.222689)
    mat.add_element("Fe",0.777311)
    return mat
end

"""
    Ferrous_Sulfate_Dosimeter_Solution()

Build the predefined NIST material **FERROUS SULFATE DOSIMETER SOLUTION** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.024` g/cm³
- Mean excitation energy ``I``: `76.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.108259`
    - `N` ⇒ `0.000027`
    - `O` ⇒ `0.878636`
    - `Na` ⇒ `0.000022`
    - `S` ⇒ `0.012968`
    - `Cl` ⇒ `0.000034`
    - `Fe` ⇒ `0.000054`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ferrous_Sulfate_Dosimeter_Solution()
    mat = Material("ferrous_sulfate_dosimeter_solution")
    mat.set_density(1.024)
    mat.set_mean_excitation_energy(76.4) # eV
    mat.add_element("H",0.108259)
    mat.add_element("N",0.000027)
    mat.add_element("O",0.878636)
    mat.add_element("Na",0.000022)
    mat.add_element("S",0.012968)
    mat.add_element("Cl",0.000034)
    mat.add_element("Fe",0.000054)
    return mat
end

"""
    Freon_12()

Build the predefined NIST material **FREON-12** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.12` g/cm³
- Mean excitation energy ``I``: `143` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.099335`
    - `F` ⇒ `0.314247`
    - `Cl` ⇒ `0.586418`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Freon_12()
    mat = Material("freon_12")
    mat.set_density(1.12)
    mat.set_mean_excitation_energy(143) # eV
    mat.add_element("C",0.099335)
    mat.add_element("F",0.314247)
    mat.add_element("Cl",0.586418)
    return mat
end

"""
    Freon_12b2()

Build the predefined NIST material **FREON-12B2** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.8` g/cm³
- Mean excitation energy ``I``: `284.9` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.057245`
    - `F` ⇒ `0.181096`
    - `Br` ⇒ `0.761659`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Freon_12b2()
    mat = Material("freon_12b2")
    mat.set_density(1.8)
    mat.set_mean_excitation_energy(284.9) # eV
    mat.add_element("C",0.057245)
    mat.add_element("F",0.181096)
    mat.add_element("Br",0.761659)
    return mat
end

"""
    Freon_13()

Build the predefined NIST material **FREON-13** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.95` g/cm³
- Mean excitation energy ``I``: `126.6` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.114983`
    - `F` ⇒ `0.545622`
    - `Cl` ⇒ `0.339396`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Freon_13()
    mat = Material("freon_13")
    mat.set_density(0.95)
    mat.set_mean_excitation_energy(126.6) # eV
    mat.add_element("C",0.114983)
    mat.add_element("F",0.545622)
    mat.add_element("Cl",0.339396)
    return mat
end

"""
    Freon_13b1()

Build the predefined NIST material **FREON-13B1** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.5` g/cm³
- Mean excitation energy ``I``: `210.5` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.080659`
    - `F` ⇒ `0.382749`
    - `Br` ⇒ `0.536592`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Freon_13b1()
    mat = Material("freon_13b1")
    mat.set_density(1.5)
    mat.set_mean_excitation_energy(210.5) # eV
    mat.add_element("C",0.080659)
    mat.add_element("F",0.382749)
    mat.add_element("Br",0.536592)
    return mat
end

"""
    Freon_13i1()

Build the predefined NIST material **FREON-13I1** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.8` g/cm³
- Mean excitation energy ``I``: `293.5` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.061309`
    - `F` ⇒ `0.290924`
    - `I` ⇒ `0.647767`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Freon_13i1()
    mat = Material("freon_13i1")
    mat.set_density(1.8)
    mat.set_mean_excitation_energy(293.5) # eV
    mat.add_element("C",0.061309)
    mat.add_element("F",0.290924)
    mat.add_element("I",0.647767)
    return mat
end

"""
    Gadolinium_Oxysulfide()

Build the predefined NIST material **GADOLINIUM OXYSULFIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.44` g/cm³
- Mean excitation energy ``I``: `493.3` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.084528`
    - `S` ⇒ `0.084690`
    - `Gd` ⇒ `0.830782`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gadolinium_Oxysulfide()
    mat = Material("gadolinium_oxysulfide")
    mat.set_density(7.44)
    mat.set_mean_excitation_energy(493.3) # eV
    mat.add_element("O",0.084528)
    mat.add_element("S",0.084690)
    mat.add_element("Gd",0.830782)
    return mat
end

"""
    Gallium_Arsenide()

Build the predefined NIST material **GALLIUM ARSENIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.31` g/cm³
- Mean excitation energy ``I``: `384.9` eV
- Composition (element ⇒ mass weight fraction):
    - `Ga` ⇒ `0.482019`
    - `As` ⇒ `0.517981`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gallium_Arsenide()
    mat = Material("gallium_arsenide")
    mat.set_density(5.31)
    mat.set_mean_excitation_energy(384.9) # eV
    mat.add_element("Ga",0.482019)
    mat.add_element("As",0.517981)
    return mat
end

"""
    Gel_In_Photographic_Emulsion()

Build the predefined NIST material **GEL IN PHOTOGRAPHIC EMULSION** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2914` g/cm³
- Mean excitation energy ``I``: `74.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.081180`
    - `C` ⇒ `0.416060`
    - `N` ⇒ `0.111240`
    - `O` ⇒ `0.380640`
    - `S` ⇒ `0.010880`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gel_In_Photographic_Emulsion()
    mat = Material("gel_in_photographic_emulsion")
    mat.set_density(1.2914)
    mat.set_mean_excitation_energy(74.8) # eV
    mat.add_element("H",0.081180)
    mat.add_element("C",0.416060)
    mat.add_element("N",0.111240)
    mat.add_element("O",0.380640)
    mat.add_element("S",0.010880)
    return mat
end

"""
    Pyrex_Glass()

Build the predefined NIST material **Pyrex Glass** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.23` g/cm³
- Mean excitation energy ``I``: `134` eV
- Composition (element ⇒ mass weight fraction):
    - `B` ⇒ `0.040064`
    - `O` ⇒ `0.539562`
    - `Na` ⇒ `0.028191`
    - `Al` ⇒ `0.011644`
    - `Si` ⇒ `0.377220`
    - `K` ⇒ `0.003321`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Pyrex_Glass()
    mat = Material("pyrex_glass")
    mat.set_density(2.23)
    mat.set_mean_excitation_energy(134) # eV
    mat.add_element("B",0.040064)
    mat.add_element("O",0.539562)
    mat.add_element("Na",0.028191)
    mat.add_element("Al",0.011644)
    mat.add_element("Si",0.377220)
    mat.add_element("K",0.003321)
    return mat
end

"""
    Glass_Lead()

Build the predefined NIST material **GLASS, LEAD** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.22` g/cm³
- Mean excitation energy ``I``: `526.4` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.156453`
    - `Si` ⇒ `0.080866`
    - `Ti` ⇒ `0.008092`
    - `As` ⇒ `0.002651`
    - `Pb` ⇒ `0.751938`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Glass_Lead()
    mat = Material("glass_lead")
    mat.set_density(6.22)
    mat.set_mean_excitation_energy(526.4) # eV
    mat.add_element("O",0.156453)
    mat.add_element("Si",0.080866)
    mat.add_element("Ti",0.008092)
    mat.add_element("As",0.002651)
    mat.add_element("Pb",0.751938)
    return mat
end

"""
    Glass_Plate()

Build the predefined NIST material **GLASS, PLATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.4` g/cm³
- Mean excitation energy ``I``: `145.4` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.459800`
    - `Na` ⇒ `0.096441`
    - `Si` ⇒ `0.336553`
    - `Ca` ⇒ `0.107205`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Glass_Plate()
    mat = Material("glass_plate")
    mat.set_density(2.4)
    mat.set_mean_excitation_energy(145.4) # eV
    mat.add_element("O",0.459800)
    mat.add_element("Na",0.096441)
    mat.add_element("Si",0.336553)
    mat.add_element("Ca",0.107205)
    return mat
end

"""
    Glucose()

Build the predefined NIST material **GLUCOSE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.54` g/cm³
- Mean excitation energy ``I``: `77.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.071204`
    - `C` ⇒ `0.363652`
    - `O` ⇒ `0.565144`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Glucose()
    mat = Material("glucose")
    mat.set_density(1.54)
    mat.set_mean_excitation_energy(77.2) # eV
    mat.add_element("H",0.071204)
    mat.add_element("C",0.363652)
    mat.add_element("O",0.565144)
    return mat
end

"""
    Glutamine()

Build the predefined NIST material **GLUTAMINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.46` g/cm³
- Mean excitation energy ``I``: `73.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.068965`
    - `C` ⇒ `0.410926`
    - `N` ⇒ `0.191681`
    - `O` ⇒ `0.328427`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Glutamine()
    mat = Material("glutamine")
    mat.set_density(1.46)
    mat.set_mean_excitation_energy(73.3) # eV
    mat.add_element("H",0.068965)
    mat.add_element("C",0.410926)
    mat.add_element("N",0.191681)
    mat.add_element("O",0.328427)
    return mat
end

"""
    Glycerol()

Build the predefined NIST material **GLYCEROL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2613` g/cm³
- Mean excitation energy ``I``: `72.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.087554`
    - `C` ⇒ `0.391262`
    - `O` ⇒ `0.521185`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Glycerol()
    mat = Material("glycerol")
    mat.set_density(1.2613)
    mat.set_mean_excitation_energy(72.6) # eV
    mat.add_element("H",0.087554)
    mat.add_element("C",0.391262)
    mat.add_element("O",0.521185)
    return mat
end

"""
    Guanine()

Build the predefined NIST material **GUANINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.58` g/cm³
- Mean excitation energy ``I``: `75` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.033346`
    - `C` ⇒ `0.397380`
    - `N` ⇒ `0.463407`
    - `O` ⇒ `0.105867`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Guanine()
    mat = Material("guanine")
    mat.set_density(1.58)
    mat.set_mean_excitation_energy(75) # eV
    mat.add_element("H",0.033346)
    mat.add_element("C",0.397380)
    mat.add_element("N",0.463407)
    mat.add_element("O",0.105867)
    return mat
end

"""
    Gypsum_Plaster_Of_Paris()

Build the predefined NIST material **GYPSUM, PLASTER OF PARIS** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.32` g/cm³
- Mean excitation energy ``I``: `129.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.023416`
    - `O` ⇒ `0.557572`
    - `S` ⇒ `0.186215`
    - `Ca` ⇒ `0.232797`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Gypsum_Plaster_Of_Paris()
    mat = Material("gypsum_plaster_of_paris")
    mat.set_density(2.32)
    mat.set_mean_excitation_energy(129.7) # eV
    mat.add_element("H",0.023416)
    mat.add_element("O",0.557572)
    mat.add_element("S",0.186215)
    mat.add_element("Ca",0.232797)
    return mat
end

"""
    N_Heptane()

Build the predefined NIST material **N-HEPTANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.68376` g/cm³
- Mean excitation energy ``I``: `54.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.160937`
    - `C` ⇒ `0.839063`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_Heptane()
    mat = Material("n_heptane")
    mat.set_density(0.68376)
    mat.set_mean_excitation_energy(54.4) # eV
    mat.add_element("H",0.160937)
    mat.add_element("C",0.839063)
    return mat
end

"""
    N_Hexane()

Build the predefined NIST material **N-HEXANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.6603` g/cm³
- Mean excitation energy ``I``: `54` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.163741`
    - `C` ⇒ `0.836259`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_Hexane()
    mat = Material("n_hexane")
    mat.set_density(0.6603)
    mat.set_mean_excitation_energy(54) # eV
    mat.add_element("H",0.163741)
    mat.add_element("C",0.836259)
    return mat
end

"""
    Kapton_Polyimide_Film()

Build the predefined NIST material **KAPTON POLYIMIDE FILM** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.42` g/cm³
- Mean excitation energy ``I``: `79.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.026362`
    - `C` ⇒ `0.691133`
    - `N` ⇒ `0.073270`
    - `O` ⇒ `0.209235`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Kapton_Polyimide_Film()
    mat = Material("kapton_polyimide_film")
    mat.set_density(1.42)
    mat.set_mean_excitation_energy(79.6) # eV
    mat.add_element("H",0.026362)
    mat.add_element("C",0.691133)
    mat.add_element("N",0.073270)
    mat.add_element("O",0.209235)
    return mat
end

"""
    Lanthanum_Oxybromide()

Build the predefined NIST material **LANTHANUM OXYBROMIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.28` g/cm³
- Mean excitation energy ``I``: `439.7` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.068138`
    - `Br` ⇒ `0.340294`
    - `La` ⇒ `0.591568`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lanthanum_Oxybromide()
    mat = Material("lanthanum_oxybromide")
    mat.set_density(6.28)
    mat.set_mean_excitation_energy(439.7) # eV
    mat.add_element("O",0.068138)
    mat.add_element("Br",0.340294)
    mat.add_element("La",0.591568)
    return mat
end

"""
    Lanthanum_Oxysulfide()

Build the predefined NIST material **LANTHANUM OXYSULFIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.86` g/cm³
- Mean excitation energy ``I``: `421.2` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.093600`
    - `S` ⇒ `0.093778`
    - `La` ⇒ `0.812622`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lanthanum_Oxysulfide()
    mat = Material("lanthanum_oxysulfide")
    mat.set_density(5.86)
    mat.set_mean_excitation_energy(421.2) # eV
    mat.add_element("O",0.093600)
    mat.add_element("S",0.093778)
    mat.add_element("La",0.812622)
    return mat
end

"""
    Lead_Oxide()

Build the predefined NIST material **LEAD OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `9.53` g/cm³
- Mean excitation energy ``I``: `766.7` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.071682`
    - `Pb` ⇒ `0.928318`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lead_Oxide()
    mat = Material("lead_oxide")
    mat.set_density(9.53)
    mat.set_mean_excitation_energy(766.7) # eV
    mat.add_element("O",0.071682)
    mat.add_element("Pb",0.928318)
    return mat
end

"""
    Lithium_Amide()

Build the predefined NIST material **LITHIUM AMIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.178` g/cm³
- Mean excitation energy ``I``: `55.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.087783`
    - `Li` ⇒ `0.302262`
    - `N` ⇒ `0.609955`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Amide()
    mat = Material("lithium_amide")
    mat.set_density(1.178)
    mat.set_mean_excitation_energy(55.5) # eV
    mat.add_element("H",0.087783)
    mat.add_element("Li",0.302262)
    mat.add_element("N",0.609955)
    return mat
end

"""
    Lithium_Carbonate()

Build the predefined NIST material **LITHIUM CARBONATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.11` g/cm³
- Mean excitation energy ``I``: `87.9` eV
- Composition (element ⇒ mass weight fraction):
    - `Li` ⇒ `0.187871`
    - `C` ⇒ `0.162550`
    - `O` ⇒ `0.649579`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Carbonate()
    mat = Material("lithium_carbonate")
    mat.set_density(2.11)
    mat.set_mean_excitation_energy(87.9) # eV
    mat.add_element("Li",0.187871)
    mat.add_element("C",0.162550)
    mat.add_element("O",0.649579)
    return mat
end

"""
    Lithium_Fluoride()

Build the predefined NIST material **LITHIUM FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.635` g/cm³
- Mean excitation energy ``I``: `94` eV
- Composition (element ⇒ mass weight fraction):
    - `Li` ⇒ `0.267585`
    - `F` ⇒ `0.732415`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Fluoride()
    mat = Material("lithium_fluoride")
    mat.set_density(2.635)
    mat.set_mean_excitation_energy(94) # eV
    mat.add_element("Li",0.267585)
    mat.add_element("F",0.732415)
    return mat
end

"""
    Lithium_Hydride()

Build the predefined NIST material **LITHIUM HYDRIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.82` g/cm³
- Mean excitation energy ``I``: `36.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.126797`
    - `Li` ⇒ `0.873203`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Hydride()
    mat = Material("lithium_hydride")
    mat.set_density(0.82)
    mat.set_mean_excitation_energy(36.5) # eV
    mat.add_element("H",0.126797)
    mat.add_element("Li",0.873203)
    return mat
end

"""
    Lithium_Iodide()

Build the predefined NIST material **LITHIUM IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.494` g/cm³
- Mean excitation energy ``I``: `485.1` eV
- Composition (element ⇒ mass weight fraction):
    - `Li` ⇒ `0.051858`
    - `I` ⇒ `0.948142`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Iodide()
    mat = Material("lithium_iodide")
    mat.set_density(3.494)
    mat.set_mean_excitation_energy(485.1) # eV
    mat.add_element("Li",0.051858)
    mat.add_element("I",0.948142)
    return mat
end

"""
    Lithium_Oxide()

Build the predefined NIST material **LITHIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.013` g/cm³
- Mean excitation energy ``I``: `73.6` eV
- Composition (element ⇒ mass weight fraction):
    - `Li` ⇒ `0.464570`
    - `O` ⇒ `0.535430`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Oxide()
    mat = Material("lithium_oxide")
    mat.set_density(2.013)
    mat.set_mean_excitation_energy(73.6) # eV
    mat.add_element("Li",0.464570)
    mat.add_element("O",0.535430)
    return mat
end

"""
    Lithium_Tetraborate()

Build the predefined NIST material **LITHIUM TETRABORATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.44` g/cm³
- Mean excitation energy ``I``: `94.6` eV
- Composition (element ⇒ mass weight fraction):
    - `Li` ⇒ `0.082085`
    - `B` ⇒ `0.255680`
    - `O` ⇒ `0.662235`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lithium_Tetraborate()
    mat = Material("lithium_tetraborate")
    mat.set_density(2.44)
    mat.set_mean_excitation_energy(94.6) # eV
    mat.add_element("Li",0.082085)
    mat.add_element("B",0.255680)
    mat.add_element("O",0.662235)
    return mat
end

"""
    Lung_Icrp()

Build the predefined NIST material **LUNG (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.05` g/cm³
- Mean excitation energy ``I``: `75.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101278`
    - `C` ⇒ `0.102310`
    - `N` ⇒ `0.028650`
    - `O` ⇒ `0.757072`
    - `Na` ⇒ `0.001840`
    - `Mg` ⇒ `0.000730`
    - `P` ⇒ `0.000800`
    - `S` ⇒ `0.002250`
    - `Cl` ⇒ `0.002660`
    - `K` ⇒ `0.001940`
    - `Ca` ⇒ `0.000090`
    - `Fe` ⇒ `0.000370`
    - `Zn` ⇒ `0.000010`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Lung_Icrp()
    mat = Material("lung_icrp")
    mat.set_density(1.05)
    mat.set_mean_excitation_energy(75.3) # eV
    mat.add_element("H",0.101278)
    mat.add_element("C",0.102310)
    mat.add_element("N",0.028650)
    mat.add_element("O",0.757072)
    mat.add_element("Na",0.001840)
    mat.add_element("Mg",0.000730)
    mat.add_element("P",0.000800)
    mat.add_element("S",0.002250)
    mat.add_element("Cl",0.002660)
    mat.add_element("K",0.001940)
    mat.add_element("Ca",0.000090)
    mat.add_element("Fe",0.000370)
    mat.add_element("Zn",0.000010)
    return mat
end

"""
    M3_Wax()

Build the predefined NIST material **M3 WAX** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.05` g/cm³
- Mean excitation energy ``I``: `67.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.114318`
    - `C` ⇒ `0.655823`
    - `O` ⇒ `0.092183`
    - `Mg` ⇒ `0.134792`
    - `Ca` ⇒ `0.002883`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function M3_Wax()
    mat = Material("m3_wax")
    mat.set_density(1.05)
    mat.set_mean_excitation_energy(67.9) # eV
    mat.add_element("H",0.114318)
    mat.add_element("C",0.655823)
    mat.add_element("O",0.092183)
    mat.add_element("Mg",0.134792)
    mat.add_element("Ca",0.002883)
    return mat
end

"""
    Magnesium_Carbonate()

Build the predefined NIST material **MAGNESIUM CARBONATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.958` g/cm³
- Mean excitation energy ``I``: `118` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.142455`
    - `O` ⇒ `0.569278`
    - `Mg` ⇒ `0.288267`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Magnesium_Carbonate()
    mat = Material("magnesium_carbonate")
    mat.set_density(2.958)
    mat.set_mean_excitation_energy(118) # eV
    mat.add_element("C",0.142455)
    mat.add_element("O",0.569278)
    mat.add_element("Mg",0.288267)
    return mat
end

"""
    Magnesium_Fluoride()

Build the predefined NIST material **MAGNESIUM FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3` g/cm³
- Mean excitation energy ``I``: `134.3` eV
- Composition (element ⇒ mass weight fraction):
    - `F` ⇒ `0.609883`
    - `Mg` ⇒ `0.390117`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Magnesium_Fluoride()
    mat = Material("magnesium_fluoride")
    mat.set_density(3)
    mat.set_mean_excitation_energy(134.3) # eV
    mat.add_element("F",0.609883)
    mat.add_element("Mg",0.390117)
    return mat
end

"""
    Magnesium_Oxide()

Build the predefined NIST material **MAGNESIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.58` g/cm³
- Mean excitation energy ``I``: `143.8` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.396964`
    - `Mg` ⇒ `0.603036`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Magnesium_Oxide()
    mat = Material("magnesium_oxide")
    mat.set_density(3.58)
    mat.set_mean_excitation_energy(143.8) # eV
    mat.add_element("O",0.396964)
    mat.add_element("Mg",0.603036)
    return mat
end

"""
    Magnesium_Tetraborate()

Build the predefined NIST material **MAGNESIUM TETRABORATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.53` g/cm³
- Mean excitation energy ``I``: `108.3` eV
- Composition (element ⇒ mass weight fraction):
    - `B` ⇒ `0.240837`
    - `O` ⇒ `0.623790`
    - `Mg` ⇒ `0.135373`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Magnesium_Tetraborate()
    mat = Material("magnesium_tetraborate")
    mat.set_density(2.53)
    mat.set_mean_excitation_energy(108.3) # eV
    mat.add_element("B",0.240837)
    mat.add_element("O",0.623790)
    mat.add_element("Mg",0.135373)
    return mat
end

"""
    Mercuric_Iodide()

Build the predefined NIST material **MERCURIC IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.36` g/cm³
- Mean excitation energy ``I``: `684.5` eV
- Composition (element ⇒ mass weight fraction):
    - `I` ⇒ `0.558560`
    - `Hg` ⇒ `0.441440`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Mercuric_Iodide()
    mat = Material("mercuric_iodide")
    mat.set_density(6.36)
    mat.set_mean_excitation_energy(684.5) # eV
    mat.add_element("I",0.558560)
    mat.add_element("Hg",0.441440)
    return mat
end

"""
    Methane()

Build the predefined NIST material **METHANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.000667151` g/cm³
- Mean excitation energy ``I``: `41.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.251306`
    - `C` ⇒ `0.748694`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Methane()
    mat = Material("methane")
    mat.set_density(0.000667151)
    mat.set_mean_excitation_energy(41.7) # eV
    mat.add_element("H",0.251306)
    mat.add_element("C",0.748694)
    return mat
end

"""
    Methanol()

Build the predefined NIST material **METHANOL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.7914` g/cm³
- Mean excitation energy ``I``: `67.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.125822`
    - `C` ⇒ `0.374852`
    - `O` ⇒ `0.499326`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Methanol()
    mat = Material("methanol")
    mat.set_density(0.7914)
    mat.set_mean_excitation_energy(67.6) # eV
    mat.add_element("H",0.125822)
    mat.add_element("C",0.374852)
    mat.add_element("O",0.499326)
    return mat
end

"""
    Mix_D_Wax()

Build the predefined NIST material **MIX D WAX** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.99` g/cm³
- Mean excitation energy ``I``: `60.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.134040`
    - `C` ⇒ `0.777960`
    - `O` ⇒ `0.035020`
    - `Mg` ⇒ `0.038594`
    - `Ti` ⇒ `0.014386`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Mix_D_Wax()
    mat = Material("mix_d_wax")
    mat.set_density(0.99)
    mat.set_mean_excitation_energy(60.9) # eV
    mat.add_element("H",0.134040)
    mat.add_element("C",0.777960)
    mat.add_element("O",0.035020)
    mat.add_element("Mg",0.038594)
    mat.add_element("Ti",0.014386)
    return mat
end

"""
    Ms20_Tissue_Substitute()

Build the predefined NIST material **MS20 TISSUE SUBSTITUTE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1` g/cm³
- Mean excitation energy ``I``: `75.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.081192`
    - `C` ⇒ `0.583442`
    - `N` ⇒ `0.017798`
    - `O` ⇒ `0.186381`
    - `Mg` ⇒ `0.130287`
    - `Cl` ⇒ `0.000900`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Ms20_Tissue_Substitute()
    mat = Material("ms20_tissue_substitute")
    mat.set_density(1)
    mat.set_mean_excitation_energy(75.1) # eV
    mat.add_element("H",0.081192)
    mat.add_element("C",0.583442)
    mat.add_element("N",0.017798)
    mat.add_element("O",0.186381)
    mat.add_element("Mg",0.130287)
    mat.add_element("Cl",0.000900)
    return mat
end

"""
    Muscle_Skeletal_Icrp()

Build the predefined NIST material **MUSCLE, SKELETAL (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.04` g/cm³
- Mean excitation energy ``I``: `75.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.100637`
    - `C` ⇒ `0.107830`
    - `N` ⇒ `0.027680`
    - `O` ⇒ `0.754773`
    - `Na` ⇒ `0.000750`
    - `Mg` ⇒ `0.000190`
    - `P` ⇒ `0.001800`
    - `S` ⇒ `0.002410`
    - `Cl` ⇒ `0.000790`
    - `K` ⇒ `0.003020`
    - `Ca` ⇒ `0.000030`
    - `Fe` ⇒ `0.000040`
    - `Zn` ⇒ `0.000050`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Muscle_Skeletal_Icrp()
    mat = Material("muscle_skeletal_icrp")
    mat.set_density(1.04)
    mat.set_mean_excitation_energy(75.3) # eV
    mat.add_element("H",0.100637)
    mat.add_element("C",0.107830)
    mat.add_element("N",0.027680)
    mat.add_element("O",0.754773)
    mat.add_element("Na",0.000750)
    mat.add_element("Mg",0.000190)
    mat.add_element("P",0.001800)
    mat.add_element("S",0.002410)
    mat.add_element("Cl",0.000790)
    mat.add_element("K",0.003020)
    mat.add_element("Ca",0.000030)
    mat.add_element("Fe",0.000040)
    mat.add_element("Zn",0.000050)
    return mat
end

"""
    Muscle_Striated_Icru()

Build the predefined NIST material **MUSCLE, STRIATED (ICRU)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.04` g/cm³
- Mean excitation energy ``I``: `74.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101997`
    - `C` ⇒ `0.123000`
    - `N` ⇒ `0.035000`
    - `O` ⇒ `0.729003`
    - `Na` ⇒ `0.000800`
    - `Mg` ⇒ `0.000200`
    - `P` ⇒ `0.002000`
    - `S` ⇒ `0.005000`
    - `K` ⇒ `0.003000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Muscle_Striated_Icru()
    mat = Material("muscle_striated_icru")
    mat.set_density(1.04)
    mat.set_mean_excitation_energy(74.7) # eV
    mat.add_element("H",0.101997)
    mat.add_element("C",0.123000)
    mat.add_element("N",0.035000)
    mat.add_element("O",0.729003)
    mat.add_element("Na",0.000800)
    mat.add_element("Mg",0.000200)
    mat.add_element("P",0.002000)
    mat.add_element("S",0.005000)
    mat.add_element("K",0.003000)
    return mat
end

"""
    Muscle_Equivalent_Liquid_With_Sucrose()

Build the predefined NIST material **MUSCLE-EQUIVALENT LIQUID, WITH SUCROSE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.11` g/cm³
- Mean excitation energy ``I``: `74.3` eV
- State of matter: `liquid`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.098234`
    - `C` ⇒ `0.156214`
    - `N` ⇒ `0.035451`
    - `O` ⇒ `0.710100`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Muscle_Equivalent_Liquid_With_Sucrose()
    mat = Material("muscle_equivalent_liquid_with_sucrose")
    mat.set_density(1.11)
    mat.set_mean_excitation_energy(74.3) # eV
    mat.set_state_of_matter("liquid")
    mat.add_element("H",0.098234)
    mat.add_element("C",0.156214)
    mat.add_element("N",0.035451)
    mat.add_element("O",0.710100)
    return mat
end

"""
    Muscle_Equivalent_Liquid_Without_Sucrose()

Build the predefined NIST material **MUSCLE-EQUIVALENT LIQUID, WITHOUT SUCROSE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.07` g/cm³
- Mean excitation energy ``I``: `74.2` eV
- State of matter: `liquid`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101969`
    - `C` ⇒ `0.120058`
    - `N` ⇒ `0.035451`
    - `O` ⇒ `0.742522`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Muscle_Equivalent_Liquid_Without_Sucrose()
    mat = Material("muscle_equivalent_liquid_without_sucrose")
    mat.set_density(1.07)
    mat.set_mean_excitation_energy(74.2) # eV
    mat.set_state_of_matter("liquid")
    mat.add_element("H",0.101969)
    mat.add_element("C",0.120058)
    mat.add_element("N",0.035451)
    mat.add_element("O",0.742522)
    return mat
end

"""
    Naphthalene()

Build the predefined NIST material **NAPHTHALENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.145` g/cm³
- Mean excitation energy ``I``: `68.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.062909`
    - `C` ⇒ `0.937091`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Naphthalene()
    mat = Material("naphthalene")
    mat.set_density(1.145)
    mat.set_mean_excitation_energy(68.4) # eV
    mat.add_element("H",0.062909)
    mat.add_element("C",0.937091)
    return mat
end

"""
    Nitrobenzene()

Build the predefined NIST material **NITROBENZENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.19867` g/cm³
- Mean excitation energy ``I``: `75.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.040935`
    - `C` ⇒ `0.585374`
    - `N` ⇒ `0.113773`
    - `O` ⇒ `0.259918`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nitrobenzene()
    mat = Material("nitrobenzene")
    mat.set_density(1.19867)
    mat.set_mean_excitation_energy(75.8) # eV
    mat.add_element("H",0.040935)
    mat.add_element("C",0.585374)
    mat.add_element("N",0.113773)
    mat.add_element("O",0.259918)
    return mat
end

"""
    Nitrous_Oxide()

Build the predefined NIST material **NITROUS OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00183094` g/cm³
- Mean excitation energy ``I``: `84.9` eV
- Composition (element ⇒ mass weight fraction):
    - `N` ⇒ `0.636483`
    - `O` ⇒ `0.363517`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nitrous_Oxide()
    mat = Material("nitrous_oxide")
    mat.set_density(0.00183094)
    mat.set_mean_excitation_energy(84.9) # eV
    mat.add_element("N",0.636483)
    mat.add_element("O",0.363517)
    return mat
end

"""
    Nylon_Du_Pont_Elvamide_8062()

Build the predefined NIST material **NYLON, DU PONT ELVAMIDE 8062** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.08` g/cm³
- Mean excitation energy ``I``: `64.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.103509`
    - `C` ⇒ `0.648415`
    - `N` ⇒ `0.099536`
    - `O` ⇒ `0.148539`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nylon_Du_Pont_Elvamide_8062()
    mat = Material("nylon_du_pont_elvamide_8062")
    mat.set_density(1.08)
    mat.set_mean_excitation_energy(64.3) # eV
    mat.add_element("H",0.103509)
    mat.add_element("C",0.648415)
    mat.add_element("N",0.099536)
    mat.add_element("O",0.148539)
    return mat
end

"""
    Nylon_Type_6_And_Type_6_6()

Build the predefined NIST material **NYLON, TYPE 6 AND TYPE 6/6** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.14` g/cm³
- Mean excitation energy ``I``: `63.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.097976`
    - `C` ⇒ `0.636856`
    - `N` ⇒ `0.123779`
    - `O` ⇒ `0.141389`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nylon_Type_6_And_Type_6_6()
    mat = Material("nylon_type_6_and_type_6_6")
    mat.set_density(1.14)
    mat.set_mean_excitation_energy(63.9) # eV
    mat.add_element("H",0.097976)
    mat.add_element("C",0.636856)
    mat.add_element("N",0.123779)
    mat.add_element("O",0.141389)
    return mat
end

"""
    Nylon_Type_6_10()

Build the predefined NIST material **NYLON, TYPE 6/10** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.14` g/cm³
- Mean excitation energy ``I``: `63.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.107062`
    - `C` ⇒ `0.680449`
    - `N` ⇒ `0.099189`
    - `O` ⇒ `0.113300`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nylon_Type_6_10()
    mat = Material("nylon_type_6_10")
    mat.set_density(1.14)
    mat.set_mean_excitation_energy(63.2) # eV
    mat.add_element("H",0.107062)
    mat.add_element("C",0.680449)
    mat.add_element("N",0.099189)
    mat.add_element("O",0.113300)
    return mat
end

"""
    Nylon_Type_11_Rilsan()

Build the predefined NIST material **NYLON, TYPE 11 (RILSAN)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.425` g/cm³
- Mean excitation energy ``I``: `61.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.115476`
    - `C` ⇒ `0.720819`
    - `N` ⇒ `0.076417`
    - `O` ⇒ `0.087289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Nylon_Type_11_Rilsan()
    mat = Material("nylon_type_11_rilsan")
    mat.set_density(1.425)
    mat.set_mean_excitation_energy(61.6) # eV
    mat.add_element("H",0.115476)
    mat.add_element("C",0.720819)
    mat.add_element("N",0.076417)
    mat.add_element("O",0.087289)
    return mat
end

"""
    Octane_Liquid()

Build the predefined NIST material **OCTANE, LIQUID** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.7026` g/cm³
- Mean excitation energy ``I``: `54.7` eV
- State of matter: `liquid`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.158821`
    - `C` ⇒ `0.841179`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Octane_Liquid()
    mat = Material("octane_liquid")
    mat.set_density(0.7026)
    mat.set_mean_excitation_energy(54.7) # eV
    mat.set_state_of_matter("liquid")
    mat.add_element("H",0.158821)
    mat.add_element("C",0.841179)
    return mat
end

"""
    Paraffin_Wax()

Build the predefined NIST material **PARAFFIN WAX** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.93` g/cm³
- Mean excitation energy ``I``: `55.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.148605`
    - `C` ⇒ `0.851395`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Paraffin_Wax()
    mat = Material("paraffin_wax")
    mat.set_density(0.93)
    mat.set_mean_excitation_energy(55.9) # eV
    mat.add_element("H",0.148605)
    mat.add_element("C",0.851395)
    return mat
end

"""
    N_Pentane()

Build the predefined NIST material **N-PENTANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.6262` g/cm³
- Mean excitation energy ``I``: `53.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.167635`
    - `C` ⇒ `0.832365`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_Pentane()
    mat = Material("n_pentane")
    mat.set_density(0.6262)
    mat.set_mean_excitation_energy(53.6) # eV
    mat.add_element("H",0.167635)
    mat.add_element("C",0.832365)
    return mat
end

"""
    Photographic_Emulsion()

Build the predefined NIST material **PHOTOGRAPHIC EMULSION** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.815` g/cm³
- Mean excitation energy ``I``: `331` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.014100`
    - `C` ⇒ `0.072261`
    - `N` ⇒ `0.019320`
    - `O` ⇒ `0.066101`
    - `S` ⇒ `0.001890`
    - `Br` ⇒ `0.349103`
    - `Ag` ⇒ `0.474105`
    - `I` ⇒ `0.003120`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Photographic_Emulsion()
    mat = Material("photographic_emulsion")
    mat.set_density(3.815)
    mat.set_mean_excitation_energy(331) # eV
    mat.add_element("H",0.014100)
    mat.add_element("C",0.072261)
    mat.add_element("N",0.019320)
    mat.add_element("O",0.066101)
    mat.add_element("S",0.001890)
    mat.add_element("Br",0.349103)
    mat.add_element("Ag",0.474105)
    mat.add_element("I",0.003120)
    return mat
end

"""
    Plastic_Scintillator_Vinyltoluene_Based()

Build the predefined NIST material **PLASTIC SCINTILLATOR (VINYLTOLUENE BASED)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.032` g/cm³
- Mean excitation energy ``I``: `64.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.085000`
    - `C` ⇒ `0.915000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Plastic_Scintillator_Vinyltoluene_Based()
    mat = Material("plastic_scintillator_vinyltoluene_based")
    mat.set_density(1.032)
    mat.set_mean_excitation_energy(64.7) # eV
    mat.add_element("H",0.085000)
    mat.add_element("C",0.915000)
    return mat
end

"""
    Plutonium_Dioxide()

Build the predefined NIST material **PLUTONIUM DIOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.46` g/cm³
- Mean excitation energy ``I``: `746.5` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.118055`
    - `Pu` ⇒ `0.881945`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Plutonium_Dioxide()
    mat = Material("plutonium_dioxide")
    mat.set_density(11.46)
    mat.set_mean_excitation_energy(746.5) # eV
    mat.add_element("O",0.118055)
    mat.add_element("Pu",0.881945)
    return mat
end

"""
    Polyacrylonitrile()

Build the predefined NIST material **POLYACRYLONITRILE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.17` g/cm³
- Mean excitation energy ``I``: `69.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.056983`
    - `C` ⇒ `0.679056`
    - `N` ⇒ `0.263962`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyacrylonitrile()
    mat = Material("polyacrylonitrile")
    mat.set_density(1.17)
    mat.set_mean_excitation_energy(69.6) # eV
    mat.add_element("H",0.056983)
    mat.add_element("C",0.679056)
    mat.add_element("N",0.263962)
    return mat
end

"""
    Polycarbonate_Makrolon_Lexan()

Build the predefined NIST material **POLYCARBONATE (MAKROLON, LEXAN)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.2` g/cm³
- Mean excitation energy ``I``: `73.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.055491`
    - `C` ⇒ `0.755751`
    - `O` ⇒ `0.188758`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polycarbonate_Makrolon_Lexan()
    mat = Material("polycarbonate_makrolon_lexan")
    mat.set_density(1.2)
    mat.set_mean_excitation_energy(73.1) # eV
    mat.add_element("H",0.055491)
    mat.add_element("C",0.755751)
    mat.add_element("O",0.188758)
    return mat
end

"""
    Polychlorostyrene()

Build the predefined NIST material **POLYCHLOROSTYRENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.3` g/cm³
- Mean excitation energy ``I``: `81.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.061869`
    - `C` ⇒ `0.696325`
    - `Cl` ⇒ `0.241806`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polychlorostyrene()
    mat = Material("polychlorostyrene")
    mat.set_density(1.3)
    mat.set_mean_excitation_energy(81.7) # eV
    mat.add_element("H",0.061869)
    mat.add_element("C",0.696325)
    mat.add_element("Cl",0.241806)
    return mat
end

"""
    Polyethylene()

Build the predefined NIST material **POLYETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.94` g/cm³
- Mean excitation energy ``I``: `57.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.143711`
    - `C` ⇒ `0.856289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyethylene()
    mat = Material("polyethylene")
    mat.set_density(0.94)
    mat.set_mean_excitation_energy(57.4) # eV
    mat.add_element("H",0.143711)
    mat.add_element("C",0.856289)
    return mat
end

"""
    Polyethylene_Terephthalate_Mylar()

Build the predefined NIST material **POLYETHYLENE TEREPHTHALATE (MYLAR)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.4` g/cm³
- Mean excitation energy ``I``: `78.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.041959`
    - `C` ⇒ `0.625017`
    - `O` ⇒ `0.333025`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyethylene_Terephthalate_Mylar()
    mat = Material("polyethylene_terephthalate_mylar")
    mat.set_density(1.4)
    mat.set_mean_excitation_energy(78.7) # eV
    mat.add_element("H",0.041959)
    mat.add_element("C",0.625017)
    mat.add_element("O",0.333025)
    return mat
end

"""
    Polymethyl_Methacralate_Lucite_Perspex_Plexiglass()

Build the predefined NIST material **POLYMETHYL METHACRALATE (LUCITE, PERSPEX, PLEXIGLASS)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.19` g/cm³
- Mean excitation energy ``I``: `74` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.080538`
    - `C` ⇒ `0.599848`
    - `O` ⇒ `0.319614`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polymethyl_Methacralate_Lucite_Perspex_Plexiglass()
    mat = Material("polymethyl_methacralate_lucite_perspex_plexiglass")
    mat.set_density(1.19)
    mat.set_mean_excitation_energy(74) # eV
    mat.add_element("H",0.080538)
    mat.add_element("C",0.599848)
    mat.add_element("O",0.319614)
    return mat
end

"""
    Polyoxymethylene()

Build the predefined NIST material **POLYOXYMETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.425` g/cm³
- Mean excitation energy ``I``: `77.4` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.067135`
    - `C` ⇒ `0.400017`
    - `O` ⇒ `0.532848`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyoxymethylene()
    mat = Material("polyoxymethylene")
    mat.set_density(1.425)
    mat.set_mean_excitation_energy(77.4) # eV
    mat.add_element("H",0.067135)
    mat.add_element("C",0.400017)
    mat.add_element("O",0.532848)
    return mat
end

"""
    Polypropylene()

Build the predefined NIST material **POLYPROPYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.9` g/cm³
- Mean excitation energy ``I``: `56.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.143711`
    - `C` ⇒ `0.856289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polypropylene()
    mat = Material("polypropylene")
    mat.set_density(0.9)
    mat.set_mean_excitation_energy(56.5) # eV
    mat.add_element("H",0.143711)
    mat.add_element("C",0.856289)
    return mat
end

"""
    Polystyrene()

Build the predefined NIST material **POLYSTYRENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.06` g/cm³
- Mean excitation energy ``I``: `68.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.077418`
    - `C` ⇒ `0.922582`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polystyrene()
    mat = Material("polystyrene")
    mat.set_density(1.06)
    mat.set_mean_excitation_energy(68.7) # eV
    mat.add_element("H",0.077418)
    mat.add_element("C",0.922582)
    return mat
end

"""
    Polytetrafluoroethylene_Teflon()

Build the predefined NIST material **POLYTETRAFLUOROETHYLENE (TEFLON)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.2` g/cm³
- Mean excitation energy ``I``: `99.1` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.240183`
    - `F` ⇒ `0.759817`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polytetrafluoroethylene_Teflon()
    mat = Material("polytetrafluoroethylene_teflon")
    mat.set_density(2.2)
    mat.set_mean_excitation_energy(99.1) # eV
    mat.add_element("C",0.240183)
    mat.add_element("F",0.759817)
    return mat
end

"""
    Polytrifluorochloroethylene()

Build the predefined NIST material **POLYTRIFLUOROCHLOROETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.1` g/cm³
- Mean excitation energy ``I``: `120.7` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.206250`
    - `F` ⇒ `0.489354`
    - `Cl` ⇒ `0.304395`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polytrifluorochloroethylene()
    mat = Material("polytrifluorochloroethylene")
    mat.set_density(2.1)
    mat.set_mean_excitation_energy(120.7) # eV
    mat.add_element("C",0.206250)
    mat.add_element("F",0.489354)
    mat.add_element("Cl",0.304395)
    return mat
end

"""
    Polyvinyl_Acetate()

Build the predefined NIST material **POLYVINYL ACETATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.19` g/cm³
- Mean excitation energy ``I``: `73.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.070245`
    - `C` ⇒ `0.558066`
    - `O` ⇒ `0.371689`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinyl_Acetate()
    mat = Material("polyvinyl_acetate")
    mat.set_density(1.19)
    mat.set_mean_excitation_energy(73.7) # eV
    mat.add_element("H",0.070245)
    mat.add_element("C",0.558066)
    mat.add_element("O",0.371689)
    return mat
end

"""
    Polyvinyl_Alcohol()

Build the predefined NIST material **POLYVINYL ALCOHOL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.3` g/cm³
- Mean excitation energy ``I``: `69.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.091517`
    - `C` ⇒ `0.545298`
    - `O` ⇒ `0.363185`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinyl_Alcohol()
    mat = Material("polyvinyl_alcohol")
    mat.set_density(1.3)
    mat.set_mean_excitation_energy(69.7) # eV
    mat.add_element("H",0.091517)
    mat.add_element("C",0.545298)
    mat.add_element("O",0.363185)
    return mat
end

"""
    Polyvinyl_Butyral()

Build the predefined NIST material **POLYVINYL BUTYRAL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.12` g/cm³
- Mean excitation energy ``I``: `67.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.092802`
    - `C` ⇒ `0.680561`
    - `O` ⇒ `0.226637`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinyl_Butyral()
    mat = Material("polyvinyl_butyral")
    mat.set_density(1.12)
    mat.set_mean_excitation_energy(67.2) # eV
    mat.add_element("H",0.092802)
    mat.add_element("C",0.680561)
    mat.add_element("O",0.226637)
    return mat
end

"""
    Polyvinyl_Chloride()

Build the predefined NIST material **POLYVINYL CHLORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.3` g/cm³
- Mean excitation energy ``I``: `108.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.048380`
    - `C` ⇒ `0.384360`
    - `Cl` ⇒ `0.567260`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinyl_Chloride()
    mat = Material("polyvinyl_chloride")
    mat.set_density(1.3)
    mat.set_mean_excitation_energy(108.2) # eV
    mat.add_element("H",0.048380)
    mat.add_element("C",0.384360)
    mat.add_element("Cl",0.567260)
    return mat
end

"""
    Polyvinylidene_Chloride_Saran()

Build the predefined NIST material **POLYVINYLIDENE CHLORIDE, SARAN** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.7` g/cm³
- Mean excitation energy ``I``: `134.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.020793`
    - `C` ⇒ `0.247793`
    - `Cl` ⇒ `0.731413`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinylidene_Chloride_Saran()
    mat = Material("polyvinylidene_chloride_saran")
    mat.set_density(1.7)
    mat.set_mean_excitation_energy(134.3) # eV
    mat.add_element("H",0.020793)
    mat.add_element("C",0.247793)
    mat.add_element("Cl",0.731413)
    return mat
end

"""
    Polyvinylidene_Fluoride()

Build the predefined NIST material **POLYVINYLIDENE FLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.76` g/cm³
- Mean excitation energy ``I``: `88.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.031480`
    - `C` ⇒ `0.375141`
    - `F` ⇒ `0.593379`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinylidene_Fluoride()
    mat = Material("polyvinylidene_fluoride")
    mat.set_density(1.76)
    mat.set_mean_excitation_energy(88.8) # eV
    mat.add_element("H",0.031480)
    mat.add_element("C",0.375141)
    mat.add_element("F",0.593379)
    return mat
end

"""
    Polyvinyl_Pyrrolidone()

Build the predefined NIST material **POLYVINYL PYRROLIDONE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.25` g/cm³
- Mean excitation energy ``I``: `67.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.081616`
    - `C` ⇒ `0.648407`
    - `N` ⇒ `0.126024`
    - `O` ⇒ `0.143953`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Polyvinyl_Pyrrolidone()
    mat = Material("polyvinyl_pyrrolidone")
    mat.set_density(1.25)
    mat.set_mean_excitation_energy(67.7) # eV
    mat.add_element("H",0.081616)
    mat.add_element("C",0.648407)
    mat.add_element("N",0.126024)
    mat.add_element("O",0.143953)
    return mat
end

"""
    Potassium_Iodide()

Build the predefined NIST material **POTASSIUM IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.13` g/cm³
- Mean excitation energy ``I``: `431.9` eV
- Composition (element ⇒ mass weight fraction):
    - `K` ⇒ `0.235528`
    - `I` ⇒ `0.764472`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Potassium_Iodide()
    mat = Material("potassium_iodide")
    mat.set_density(3.13)
    mat.set_mean_excitation_energy(431.9) # eV
    mat.add_element("K",0.235528)
    mat.add_element("I",0.764472)
    return mat
end

"""
    Potassium_Oxide()

Build the predefined NIST material **POTASSIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.32` g/cm³
- Mean excitation energy ``I``: `189.9` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.169852`
    - `K` ⇒ `0.830148`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Potassium_Oxide()
    mat = Material("potassium_oxide")
    mat.set_density(2.32)
    mat.set_mean_excitation_energy(189.9) # eV
    mat.add_element("O",0.169852)
    mat.add_element("K",0.830148)
    return mat
end

"""
    Propane()

Build the predefined NIST material **PROPANE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00187939` g/cm³
- Mean excitation energy ``I``: `47.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.182855`
    - `C` ⇒ `0.817145`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Propane()
    mat = Material("propane")
    mat.set_density(0.00187939)
    mat.set_mean_excitation_energy(47.1) # eV
    mat.add_element("H",0.182855)
    mat.add_element("C",0.817145)
    return mat
end

"""
    Propane_Liquid()

Build the predefined NIST material **PROPANE, LIQUID** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.43` g/cm³
- Mean excitation energy ``I``: `52` eV
- State of matter: `liquid`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.182855`
    - `C` ⇒ `0.817145`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Propane_Liquid()
    mat = Material("propane_liquid")
    mat.set_density(0.43)
    mat.set_mean_excitation_energy(52) # eV
    mat.set_state_of_matter("liquid")
    mat.add_element("H",0.182855)
    mat.add_element("C",0.817145)
    return mat
end

"""
    N_Propyl_Alcohol()

Build the predefined NIST material **N-PROPYL ALCOHOL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.8035` g/cm³
- Mean excitation energy ``I``: `61.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.134173`
    - `C` ⇒ `0.599595`
    - `O` ⇒ `0.266232`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function N_Propyl_Alcohol()
    mat = Material("n_propyl_alcohol")
    mat.set_density(0.8035)
    mat.set_mean_excitation_energy(61.1) # eV
    mat.add_element("H",0.134173)
    mat.add_element("C",0.599595)
    mat.add_element("O",0.266232)
    return mat
end

"""
    Pyridine()

Build the predefined NIST material **PYRIDINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.9819` g/cm³
- Mean excitation energy ``I``: `66.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.063710`
    - `C` ⇒ `0.759217`
    - `N` ⇒ `0.177073`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Pyridine()
    mat = Material("pyridine")
    mat.set_density(0.9819)
    mat.set_mean_excitation_energy(66.2) # eV
    mat.add_element("H",0.063710)
    mat.add_element("C",0.759217)
    mat.add_element("N",0.177073)
    return mat
end

"""
    Rubber_Butyl()

Build the predefined NIST material **RUBBER, BUTYL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.92` g/cm³
- Mean excitation energy ``I``: `56.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.143711`
    - `C` ⇒ `0.856289`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rubber_Butyl()
    mat = Material("rubber_butyl")
    mat.set_density(0.92)
    mat.set_mean_excitation_energy(56.5) # eV
    mat.add_element("H",0.143711)
    mat.add_element("C",0.856289)
    return mat
end

"""
    Rubber_Natural()

Build the predefined NIST material **RUBBER, NATURAL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.92` g/cm³
- Mean excitation energy ``I``: `59.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.118371`
    - `C` ⇒ `0.881629`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rubber_Natural()
    mat = Material("rubber_natural")
    mat.set_density(0.92)
    mat.set_mean_excitation_energy(59.8) # eV
    mat.add_element("H",0.118371)
    mat.add_element("C",0.881629)
    return mat
end

"""
    Rubber_Neoprene()

Build the predefined NIST material **RUBBER, NEOPRENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.23` g/cm³
- Mean excitation energy ``I``: `93` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.056920`
    - `C` ⇒ `0.542646`
    - `Cl` ⇒ `0.400434`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Rubber_Neoprene()
    mat = Material("rubber_neoprene")
    mat.set_density(1.23)
    mat.set_mean_excitation_energy(93) # eV
    mat.add_element("H",0.056920)
    mat.add_element("C",0.542646)
    mat.add_element("Cl",0.400434)
    return mat
end

"""
    Silicon_Dioxide()

Build the predefined NIST material **SILICON DIOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.32` g/cm³
- Mean excitation energy ``I``: `139.2` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.532565`
    - `Si` ⇒ `0.467435`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silicon_Dioxide()
    mat = Material("silicon_dioxide")
    mat.set_density(2.32)
    mat.set_mean_excitation_energy(139.2) # eV
    mat.add_element("O",0.532565)
    mat.add_element("Si",0.467435)
    return mat
end

"""
    Silver_Bromide()

Build the predefined NIST material **SILVER BROMIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.473` g/cm³
- Mean excitation energy ``I``: `486.6` eV
- Composition (element ⇒ mass weight fraction):
    - `Br` ⇒ `0.425537`
    - `Ag` ⇒ `0.574463`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silver_Bromide()
    mat = Material("silver_bromide")
    mat.set_density(6.473)
    mat.set_mean_excitation_energy(486.6) # eV
    mat.add_element("Br",0.425537)
    mat.add_element("Ag",0.574463)
    return mat
end

"""
    Silver_Chloride()

Build the predefined NIST material **SILVER CHLORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `5.56` g/cm³
- Mean excitation energy ``I``: `398.4` eV
- Composition (element ⇒ mass weight fraction):
    - `Cl` ⇒ `0.247368`
    - `Ag` ⇒ `0.752632`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silver_Chloride()
    mat = Material("silver_chloride")
    mat.set_density(5.56)
    mat.set_mean_excitation_energy(398.4) # eV
    mat.add_element("Cl",0.247368)
    mat.add_element("Ag",0.752632)
    return mat
end

"""
    Silver_Halides_In_Photographic_Emulsion()

Build the predefined NIST material **SILVER HALIDES IN PHOTOGRAPHIC EMULSION** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.47` g/cm³
- Mean excitation energy ``I``: `487.1` eV
- Composition (element ⇒ mass weight fraction):
    - `Br` ⇒ `0.422895`
    - `Ag` ⇒ `0.573748`
    - `I` ⇒ `0.003357`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silver_Halides_In_Photographic_Emulsion()
    mat = Material("silver_halides_in_photographic_emulsion")
    mat.set_density(6.47)
    mat.set_mean_excitation_energy(487.1) # eV
    mat.add_element("Br",0.422895)
    mat.add_element("Ag",0.573748)
    mat.add_element("I",0.003357)
    return mat
end

"""
    Silver_Iodide()

Build the predefined NIST material **SILVER IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `6.01` g/cm³
- Mean excitation energy ``I``: `543.5` eV
- Composition (element ⇒ mass weight fraction):
    - `Ag` ⇒ `0.459458`
    - `I` ⇒ `0.540542`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Silver_Iodide()
    mat = Material("silver_iodide")
    mat.set_density(6.01)
    mat.set_mean_excitation_energy(543.5) # eV
    mat.add_element("Ag",0.459458)
    mat.add_element("I",0.540542)
    return mat
end

"""
    Skin_Icrp()

Build the predefined NIST material **SKIN (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.1` g/cm³
- Mean excitation energy ``I``: `72.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.100588`
    - `C` ⇒ `0.228250`
    - `N` ⇒ `0.046420`
    - `O` ⇒ `0.619002`
    - `Na` ⇒ `0.000070`
    - `Mg` ⇒ `0.000060`
    - `P` ⇒ `0.000330`
    - `S` ⇒ `0.001590`
    - `Cl` ⇒ `0.002670`
    - `K` ⇒ `0.000850`
    - `Ca` ⇒ `0.000150`
    - `Fe` ⇒ `0.000010`
    - `Zn` ⇒ `0.000010`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Skin_Icrp()
    mat = Material("skin_icrp")
    mat.set_density(1.1)
    mat.set_mean_excitation_energy(72.7) # eV
    mat.add_element("H",0.100588)
    mat.add_element("C",0.228250)
    mat.add_element("N",0.046420)
    mat.add_element("O",0.619002)
    mat.add_element("Na",0.000070)
    mat.add_element("Mg",0.000060)
    mat.add_element("P",0.000330)
    mat.add_element("S",0.001590)
    mat.add_element("Cl",0.002670)
    mat.add_element("K",0.000850)
    mat.add_element("Ca",0.000150)
    mat.add_element("Fe",0.000010)
    mat.add_element("Zn",0.000010)
    return mat
end

"""
    Sodium_Carbonate()

Build the predefined NIST material **SODIUM CARBONATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.532` g/cm³
- Mean excitation energy ``I``: `125` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.113323`
    - `O` ⇒ `0.452861`
    - `Na` ⇒ `0.433815`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sodium_Carbonate()
    mat = Material("sodium_carbonate")
    mat.set_density(2.532)
    mat.set_mean_excitation_energy(125) # eV
    mat.add_element("C",0.113323)
    mat.add_element("O",0.452861)
    mat.add_element("Na",0.433815)
    return mat
end

"""
    Sodium_Iodide()

Build the predefined NIST material **SODIUM IODIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `3.667` g/cm³
- Mean excitation energy ``I``: `452` eV
- Composition (element ⇒ mass weight fraction):
    - `Na` ⇒ `0.153373`
    - `I` ⇒ `0.846627`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sodium_Iodide()
    mat = Material("sodium_iodide")
    mat.set_density(3.667)
    mat.set_mean_excitation_energy(452) # eV
    mat.add_element("Na",0.153373)
    mat.add_element("I",0.846627)
    return mat
end

"""
    Sodium_Monoxide()

Build the predefined NIST material **SODIUM MONOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.27` g/cm³
- Mean excitation energy ``I``: `148.8` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.258143`
    - `Na` ⇒ `0.741857`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sodium_Monoxide()
    mat = Material("sodium_monoxide")
    mat.set_density(2.27)
    mat.set_mean_excitation_energy(148.8) # eV
    mat.add_element("O",0.258143)
    mat.add_element("Na",0.741857)
    return mat
end

"""
    Sodium_Nitrate()

Build the predefined NIST material **SODIUM NITRATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.261` g/cm³
- Mean excitation energy ``I``: `114.6` eV
- Composition (element ⇒ mass weight fraction):
    - `N` ⇒ `0.164795`
    - `O` ⇒ `0.564720`
    - `Na` ⇒ `0.270485`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sodium_Nitrate()
    mat = Material("sodium_nitrate")
    mat.set_density(2.261)
    mat.set_mean_excitation_energy(114.6) # eV
    mat.add_element("N",0.164795)
    mat.add_element("O",0.564720)
    mat.add_element("Na",0.270485)
    return mat
end

"""
    Stilbene()

Build the predefined NIST material **STILBENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.9707` g/cm³
- Mean excitation energy ``I``: `67.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.067101`
    - `C` ⇒ `0.932899`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Stilbene()
    mat = Material("stilbene")
    mat.set_density(0.9707)
    mat.set_mean_excitation_energy(67.7) # eV
    mat.add_element("H",0.067101)
    mat.add_element("C",0.932899)
    return mat
end

"""
    Sucrose()

Build the predefined NIST material **SUCROSE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.5805` g/cm³
- Mean excitation energy ``I``: `77.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.064779`
    - `C` ⇒ `0.421070`
    - `O` ⇒ `0.514151`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Sucrose()
    mat = Material("sucrose")
    mat.set_density(1.5805)
    mat.set_mean_excitation_energy(77.5) # eV
    mat.add_element("H",0.064779)
    mat.add_element("C",0.421070)
    mat.add_element("O",0.514151)
    return mat
end

"""
    Terphenyl()

Build the predefined NIST material **TERPHENYL** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.234` g/cm³
- Mean excitation energy ``I``: `71.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.044543`
    - `C` ⇒ `0.955457`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Terphenyl()
    mat = Material("terphenyl")
    mat.set_density(1.234)
    mat.set_mean_excitation_energy(71.7) # eV
    mat.add_element("H",0.044543)
    mat.add_element("C",0.955457)
    return mat
end

"""
    Testes_Icrp()

Build the predefined NIST material **TESTES (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.04` g/cm³
- Mean excitation energy ``I``: `75` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.104166`
    - `C` ⇒ `0.092270`
    - `N` ⇒ `0.019940`
    - `O` ⇒ `0.773884`
    - `Na` ⇒ `0.002260`
    - `Mg` ⇒ `0.000110`
    - `P` ⇒ `0.001250`
    - `S` ⇒ `0.001460`
    - `Cl` ⇒ `0.002440`
    - `K` ⇒ `0.002080`
    - `Ca` ⇒ `0.000100`
    - `Fe` ⇒ `0.000020`
    - `Zn` ⇒ `0.000020`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Testes_Icrp()
    mat = Material("testes_icrp")
    mat.set_density(1.04)
    mat.set_mean_excitation_energy(75) # eV
    mat.add_element("H",0.104166)
    mat.add_element("C",0.092270)
    mat.add_element("N",0.019940)
    mat.add_element("O",0.773884)
    mat.add_element("Na",0.002260)
    mat.add_element("Mg",0.000110)
    mat.add_element("P",0.001250)
    mat.add_element("S",0.001460)
    mat.add_element("Cl",0.002440)
    mat.add_element("K",0.002080)
    mat.add_element("Ca",0.000100)
    mat.add_element("Fe",0.000020)
    mat.add_element("Zn",0.000020)
    return mat
end

"""
    Tetrachloroethylene()

Build the predefined NIST material **TETRACHLOROETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.625` g/cm³
- Mean excitation energy ``I``: `159.2` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.144856`
    - `Cl` ⇒ `0.855144`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tetrachloroethylene()
    mat = Material("tetrachloroethylene")
    mat.set_density(1.625)
    mat.set_mean_excitation_energy(159.2) # eV
    mat.add_element("C",0.144856)
    mat.add_element("Cl",0.855144)
    return mat
end

"""
    Thallium_Chloride()

Build the predefined NIST material **THALLIUM CHLORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `7.004` g/cm³
- Mean excitation energy ``I``: `690.3` eV
- Composition (element ⇒ mass weight fraction):
    - `Cl` ⇒ `0.147822`
    - `Tl` ⇒ `0.852178`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Thallium_Chloride()
    mat = Material("thallium_chloride")
    mat.set_density(7.004)
    mat.set_mean_excitation_energy(690.3) # eV
    mat.add_element("Cl",0.147822)
    mat.add_element("Tl",0.852178)
    return mat
end

"""
    Tissue_Soft_Icrp()

Build the predefined NIST material **TISSUE, SOFT (ICRP)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1` g/cm³
- Mean excitation energy ``I``: `72.3` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.104472`
    - `C` ⇒ `0.232190`
    - `N` ⇒ `0.024880`
    - `O` ⇒ `0.630238`
    - `Na` ⇒ `0.001130`
    - `Mg` ⇒ `0.000130`
    - `P` ⇒ `0.001330`
    - `S` ⇒ `0.001990`
    - `Cl` ⇒ `0.001340`
    - `K` ⇒ `0.001990`
    - `Ca` ⇒ `0.000230`
    - `Fe` ⇒ `0.000050`
    - `Zn` ⇒ `0.000030`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tissue_Soft_Icrp()
    mat = Material("tissue_soft_icrp")
    mat.set_density(1)
    mat.set_mean_excitation_energy(72.3) # eV
    mat.add_element("H",0.104472)
    mat.add_element("C",0.232190)
    mat.add_element("N",0.024880)
    mat.add_element("O",0.630238)
    mat.add_element("Na",0.001130)
    mat.add_element("Mg",0.000130)
    mat.add_element("P",0.001330)
    mat.add_element("S",0.001990)
    mat.add_element("Cl",0.001340)
    mat.add_element("K",0.001990)
    mat.add_element("Ca",0.000230)
    mat.add_element("Fe",0.000050)
    mat.add_element("Zn",0.000030)
    return mat
end

"""
    Tissue_Soft_Icru_Four_Component()

Build the predefined NIST material **TISSUE, SOFT (ICRU FOUR-COMPONENT)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1` g/cm³
- Mean excitation energy ``I``: `74.9` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101172`
    - `C` ⇒ `0.111000`
    - `N` ⇒ `0.026000`
    - `O` ⇒ `0.761828`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tissue_Soft_Icru_Four_Component()
    mat = Material("tissue_soft_icru_four_component")
    mat.set_density(1)
    mat.set_mean_excitation_energy(74.9) # eV
    mat.add_element("H",0.101172)
    mat.add_element("C",0.111000)
    mat.add_element("N",0.026000)
    mat.add_element("O",0.761828)
    return mat
end

"""
    Tissue_Equivalent_Gas_Methane_Based()

Build the predefined NIST material **TISSUE-EQUIVALENT GAS (METHANE BASED)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00106409` g/cm³
- Mean excitation energy ``I``: `61.2` eV
- State of matter: `gaz`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.101869`
    - `C` ⇒ `0.456179`
    - `N` ⇒ `0.035172`
    - `O` ⇒ `0.406780`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tissue_Equivalent_Gas_Methane_Based()
    mat = Material("tissue_equivalent_gas_methane_based")
    mat.set_density(0.00106409)
    mat.set_mean_excitation_energy(61.2) # eV
    mat.set_state_of_matter("gaz")
    mat.add_element("H",0.101869)
    mat.add_element("C",0.456179)
    mat.add_element("N",0.035172)
    mat.add_element("O",0.406780)
    return mat
end

"""
    Tissue_Equivalent_Gas_Propane_Based()

Build the predefined NIST material **TISSUE-EQUIVALENT GAS (PROPANE BASED)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.00182628` g/cm³
- Mean excitation energy ``I``: `59.5` eV
- State of matter: `gaz`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.102672`
    - `C` ⇒ `0.568940`
    - `N` ⇒ `0.035022`
    - `O` ⇒ `0.293366`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tissue_Equivalent_Gas_Propane_Based()
    mat = Material("tissue_equivalent_gas_propane_based")
    mat.set_density(0.00182628)
    mat.set_mean_excitation_energy(59.5) # eV
    mat.set_state_of_matter("gaz")
    mat.add_element("H",0.102672)
    mat.add_element("C",0.568940)
    mat.add_element("N",0.035022)
    mat.add_element("O",0.293366)
    return mat
end

"""
    Titanium_Dioxide()

Build the predefined NIST material **TITANIUM DIOXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `4.26` g/cm³
- Mean excitation energy ``I``: `179.5` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.400592`
    - `Ti` ⇒ `0.599408`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Titanium_Dioxide()
    mat = Material("titanium_dioxide")
    mat.set_density(4.26)
    mat.set_mean_excitation_energy(179.5) # eV
    mat.add_element("O",0.400592)
    mat.add_element("Ti",0.599408)
    return mat
end

"""
    Toluene()

Build the predefined NIST material **TOLUENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.8669` g/cm³
- Mean excitation energy ``I``: `62.5` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.087510`
    - `C` ⇒ `0.912490`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Toluene()
    mat = Material("toluene")
    mat.set_density(0.8669)
    mat.set_mean_excitation_energy(62.5) # eV
    mat.add_element("H",0.087510)
    mat.add_element("C",0.912490)
    return mat
end

"""
    Trichloroethylene()

Build the predefined NIST material **TRICHLOROETHYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.46` g/cm³
- Mean excitation energy ``I``: `148.1` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.007671`
    - `C` ⇒ `0.182831`
    - `Cl` ⇒ `0.809498`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Trichloroethylene()
    mat = Material("trichloroethylene")
    mat.set_density(1.46)
    mat.set_mean_excitation_energy(148.1) # eV
    mat.add_element("H",0.007671)
    mat.add_element("C",0.182831)
    mat.add_element("Cl",0.809498)
    return mat
end

"""
    Triethyl_Phosphate()

Build the predefined NIST material **TRIETHYL PHOSPHATE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.07` g/cm³
- Mean excitation energy ``I``: `81.2` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.082998`
    - `C` ⇒ `0.395628`
    - `O` ⇒ `0.351334`
    - `P` ⇒ `0.170040`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Triethyl_Phosphate()
    mat = Material("triethyl_phosphate")
    mat.set_density(1.07)
    mat.set_mean_excitation_energy(81.2) # eV
    mat.add_element("H",0.082998)
    mat.add_element("C",0.395628)
    mat.add_element("O",0.351334)
    mat.add_element("P",0.170040)
    return mat
end

"""
    Tungsten_Hexafluoride()

Build the predefined NIST material **TUNGSTEN HEXAFLUORIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `2.4` g/cm³
- Mean excitation energy ``I``: `354.4` eV
- Composition (element ⇒ mass weight fraction):
    - `F` ⇒ `0.382723`
    - `W` ⇒ `0.617277`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Tungsten_Hexafluoride()
    mat = Material("tungsten_hexafluoride")
    mat.set_density(2.4)
    mat.set_mean_excitation_energy(354.4) # eV
    mat.add_element("F",0.382723)
    mat.add_element("W",0.617277)
    return mat
end

"""
    Uranium_Dicarbide()

Build the predefined NIST material **URANIUM DICARBIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `11.28` g/cm³
- Mean excitation energy ``I``: `752` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.091669`
    - `U` ⇒ `0.908331`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Uranium_Dicarbide()
    mat = Material("uranium_dicarbide")
    mat.set_density(11.28)
    mat.set_mean_excitation_energy(752) # eV
    mat.add_element("C",0.091669)
    mat.add_element("U",0.908331)
    return mat
end

"""
    Uranium_Monocarbide()

Build the predefined NIST material **URANIUM MONOCARBIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `13.63` g/cm³
- Mean excitation energy ``I``: `862` eV
- Composition (element ⇒ mass weight fraction):
    - `C` ⇒ `0.048036`
    - `U` ⇒ `0.951964`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Uranium_Monocarbide()
    mat = Material("uranium_monocarbide")
    mat.set_density(13.63)
    mat.set_mean_excitation_energy(862) # eV
    mat.add_element("C",0.048036)
    mat.add_element("U",0.951964)
    return mat
end

"""
    Uranium_Oxide()

Build the predefined NIST material **URANIUM OXIDE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `10.96` g/cm³
- Mean excitation energy ``I``: `720.6` eV
- Composition (element ⇒ mass weight fraction):
    - `O` ⇒ `0.118502`
    - `U` ⇒ `0.881498`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Uranium_Oxide()
    mat = Material("uranium_oxide")
    mat.set_density(10.96)
    mat.set_mean_excitation_energy(720.6) # eV
    mat.add_element("O",0.118502)
    mat.add_element("U",0.881498)
    return mat
end

"""
    Urea()

Build the predefined NIST material **UREA** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.323` g/cm³
- Mean excitation energy ``I``: `72.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.067131`
    - `C` ⇒ `0.199999`
    - `N` ⇒ `0.466459`
    - `O` ⇒ `0.266411`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Urea()
    mat = Material("urea")
    mat.set_density(1.323)
    mat.set_mean_excitation_energy(72.8) # eV
    mat.add_element("H",0.067131)
    mat.add_element("C",0.199999)
    mat.add_element("N",0.466459)
    mat.add_element("O",0.266411)
    return mat
end

"""
    Valine()

Build the predefined NIST material **VALINE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.23` g/cm³
- Mean excitation energy ``I``: `67.7` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.094641`
    - `C` ⇒ `0.512645`
    - `N` ⇒ `0.119565`
    - `O` ⇒ `0.273150`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Valine()
    mat = Material("valine")
    mat.set_density(1.23)
    mat.set_mean_excitation_energy(67.7) # eV
    mat.add_element("H",0.094641)
    mat.add_element("C",0.512645)
    mat.add_element("N",0.119565)
    mat.add_element("O",0.273150)
    return mat
end

"""
    Viton_Fluoroelastomer()

Build the predefined NIST material **VITON FLUOROELASTOMER** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.8` g/cm³
- Mean excitation energy ``I``: `98.6` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.009417`
    - `C` ⇒ `0.280555`
    - `F` ⇒ `0.710028`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Viton_Fluoroelastomer()
    mat = Material("viton_fluoroelastomer")
    mat.set_density(1.8)
    mat.set_mean_excitation_energy(98.6) # eV
    mat.add_element("H",0.009417)
    mat.add_element("C",0.280555)
    mat.add_element("F",0.710028)
    return mat
end

"""
    Water()

Build the predefined NIST material **WATER, LIQUID** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1` g/cm³
- Mean excitation energy ``I``: `75` eV
- State of matter: `liquid`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.111894`
    - `O` ⇒ `0.888106`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Water()
    mat = Material("water")
    mat.set_density(1)
    mat.set_mean_excitation_energy(75) # eV
    mat.set_state_of_matter("liquid")
    mat.add_element("H",0.111894)
    mat.add_element("O",0.888106)
    return mat
end

"""
    Water_Vapor()

Build the predefined NIST material **WATER VAPOR** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.000756182` g/cm³
- Mean excitation energy ``I``: `71.6` eV
- State of matter: `gaz`
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.111894`
    - `O` ⇒ `0.888106`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Water_Vapor()
    mat = Material("water_vapor")
    mat.set_density(0.000756182)
    mat.set_mean_excitation_energy(71.6) # eV
    mat.set_state_of_matter("gaz")
    mat.add_element("H",0.111894)
    mat.add_element("O",0.888106)
    return mat
end

"""
    Xylene()

Build the predefined NIST material **XYLENE** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `0.87` g/cm³
- Mean excitation energy ``I``: `61.8` eV
- Composition (element ⇒ mass weight fraction):
    - `H` ⇒ `0.094935`
    - `C` ⇒ `0.905065`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Xylene()
    mat = Material("xylene")
    mat.set_density(0.87)
    mat.set_mean_excitation_energy(61.8) # eV
    mat.add_element("H",0.094935)
    mat.add_element("C",0.905065)
    return mat
end

"""
    Carbon_Graphite()

Build the predefined NIST material **CARBON (GRAPHITE)** and return it as a ready-to-use
`Material` object with the following preset properties:

- Mass density: `1.7` g/cm³
- Mean excitation energy ``I``: `78` eV
- Composition (element ⇒ mass weight fraction): `C` ⇒ `1.000000`

Any property can be overridden after construction (for instance with
`set_density`, `set_mean_excitation_energy`, `set_state_of_matter` or
`add_element`).
"""
function Carbon_Graphite()
    mat = Material("carbon_graphite")
    mat.set_density(1.7)
    mat.set_mean_excitation_energy(78) # eV
    mat.add_element("C",1.000000)
    return mat
end

const MATERIAL_CONSTRUCTORS = [
    :Hydrogen,
    :Helium,
    :Lithium,
    :Beryllium,
    :Boron,
    :Amorphous_Carbon,
    :Nitrogen,
    :Oxygen,
    :Fluorine,
    :Neon,
    :Sodium,
    :Magnesium,
    :Aluminum,
    :Silicon,
    :Phosphorus,
    :Sulfur,
    :Chlorine,
    :Argon,
    :Potassium,
    :Calcium,
    :Scandium,
    :Titanium,
    :Vanadium,
    :Chromium,
    :Manganese,
    :Iron,
    :Cobalt,
    :Nickel,
    :Copper,
    :Zinc,
    :Gallium,
    :Germanium,
    :Arsenic,
    :Selenium,
    :Bromine,
    :Krypton,
    :Rubidium,
    :Strontium,
    :Yttrium,
    :Zirconium,
    :Niobium,
    :Molybdenum,
    :Technetium,
    :Ruthenium,
    :Rhodium,
    :Palladium,
    :Silver,
    :Cadmium,
    :Indium,
    :Tin,
    :Antimony,
    :Tellurium,
    :Iodine,
    :Xenon,
    :Cesium,
    :Barium,
    :Lanthanum,
    :Cerium,
    :Praseodymium,
    :Neodymium,
    :Promethium,
    :Samarium,
    :Europium,
    :Gadolinium,
    :Terbium,
    :Dysprosium,
    :Holmium,
    :Erbium,
    :Thulium,
    :Ytterbium,
    :Lutetium,
    :Hafnium,
    :Tantalum,
    :Tungsten,
    :Rhenium,
    :Osmium,
    :Iridium,
    :Platinum,
    :Gold,
    :Mercury,
    :Thallium,
    :Lead,
    :Bismuth,
    :Polonium,
    :Astatine,
    :Radon,
    :Francium,
    :Radium,
    :Actinium,
    :Thorium,
    :Protactinium,
    :Uranium,
    :Neptunium,
    :Plutonium,
    :Americium,
    :Curium,
    :Berkelium,
    :Californium,
    :A_150_Tissue_Equivalent_Plastic,
    :Acetone,
    :Acetylene,
    :Adenine,
    :Adipose_Tissue_Icrp,
    :Air_Dry_Near_Sea_Level,
    :Alanine,
    :Aluminum_Oxide,
    :Amber,
    :Ammonia,
    :Aniline,
    :Anthracene,
    :B_100_Bone_Equivalent_Plastic,
    :Bakelite,
    :Barium_Fluoride,
    :Barium_Sulfate,
    :Benzene,
    :Beryllium_Oxide,
    :Bismuth_Germanium_Oxide,
    :Blood_Icrp,
    :Bone_Compact_Icru,
    :Bone_Cortical_Icrp,
    :Boron_Carbide,
    :Boron_Oxide,
    :Brain_Icrp,
    :Butane,
    :N_Butyl_Alcohol,
    :C_552_Air_Equivalent_Plastic,
    :Cadmium_Telluride,
    :Cadmium_Tungstate,
    :Calcium_Carbonate,
    :Calcium_Fluoride,
    :Calcium_Oxide,
    :Calcium_Sulfate,
    :Calcium_Tungstate,
    :Carbon_Dioxide,
    :Carbon_Tetrachloride,
    :Cellulose_Acetate_Cellophane,
    :Cellulose_Acetate_Butyrate,
    :Cellulose_Nitrate,
    :Ceric_Sulfate_Dosimeter_Solution,
    :Cesium_Fluoride,
    :Cesium_Iodide,
    :Chlorobenzene,
    :Chloroform,
    :Concrete_Portland,
    :Cyclohexane,
    :M1_2_Dichlorobenzene,
    :Dichlorodiethyl_Ether,
    :M1_2_Dichloroethane,
    :Diethyl_Ether,
    :N_N_Dimethyl_Formamide,
    :Dimethyl_Sulfoxide,
    :Ethane,
    :Ethyl_Alcohol,
    :Ethyl_Cellulose,
    :Ethylene,
    :Eye_Lens_Icrp,
    :Ferric_Oxide,
    :Ferroboride,
    :Ferrous_Oxide,
    :Ferrous_Sulfate_Dosimeter_Solution,
    :Freon_12,
    :Freon_12b2,
    :Freon_13,
    :Freon_13b1,
    :Freon_13i1,
    :Gadolinium_Oxysulfide,
    :Gallium_Arsenide,
    :Gel_In_Photographic_Emulsion,
    :Pyrex_Glass,
    :Glass_Lead,
    :Glass_Plate,
    :Glucose,
    :Glutamine,
    :Glycerol,
    :Guanine,
    :Gypsum_Plaster_Of_Paris,
    :N_Heptane,
    :N_Hexane,
    :Kapton_Polyimide_Film,
    :Lanthanum_Oxybromide,
    :Lanthanum_Oxysulfide,
    :Lead_Oxide,
    :Lithium_Amide,
    :Lithium_Carbonate,
    :Lithium_Fluoride,
    :Lithium_Hydride,
    :Lithium_Iodide,
    :Lithium_Oxide,
    :Lithium_Tetraborate,
    :Lung_Icrp,
    :M3_Wax,
    :Magnesium_Carbonate,
    :Magnesium_Fluoride,
    :Magnesium_Oxide,
    :Magnesium_Tetraborate,
    :Mercuric_Iodide,
    :Methane,
    :Methanol,
    :Mix_D_Wax,
    :Ms20_Tissue_Substitute,
    :Muscle_Skeletal_Icrp,
    :Muscle_Striated_Icru,
    :Muscle_Equivalent_Liquid_With_Sucrose,
    :Muscle_Equivalent_Liquid_Without_Sucrose,
    :Naphthalene,
    :Nitrobenzene,
    :Nitrous_Oxide,
    :Nylon_Du_Pont_Elvamide_8062,
    :Nylon_Type_6_And_Type_6_6,
    :Nylon_Type_6_10,
    :Nylon_Type_11_Rilsan,
    :Octane_Liquid,
    :Paraffin_Wax,
    :N_Pentane,
    :Photographic_Emulsion,
    :Plastic_Scintillator_Vinyltoluene_Based,
    :Plutonium_Dioxide,
    :Polyacrylonitrile,
    :Polycarbonate_Makrolon_Lexan,
    :Polychlorostyrene,
    :Polyethylene,
    :Polyethylene_Terephthalate_Mylar,
    :Polymethyl_Methacralate_Lucite_Perspex_Plexiglass,
    :Polyoxymethylene,
    :Polypropylene,
    :Polystyrene,
    :Polytetrafluoroethylene_Teflon,
    :Polytrifluorochloroethylene,
    :Polyvinyl_Acetate,
    :Polyvinyl_Alcohol,
    :Polyvinyl_Butyral,
    :Polyvinyl_Chloride,
    :Polyvinylidene_Chloride_Saran,
    :Polyvinylidene_Fluoride,
    :Polyvinyl_Pyrrolidone,
    :Potassium_Iodide,
    :Potassium_Oxide,
    :Propane,
    :Propane_Liquid,
    :N_Propyl_Alcohol,
    :Pyridine,
    :Rubber_Butyl,
    :Rubber_Natural,
    :Rubber_Neoprene,
    :Silicon_Dioxide,
    :Silver_Bromide,
    :Silver_Chloride,
    :Silver_Halides_In_Photographic_Emulsion,
    :Silver_Iodide,
    :Skin_Icrp,
    :Sodium_Carbonate,
    :Sodium_Iodide,
    :Sodium_Monoxide,
    :Sodium_Nitrate,
    :Stilbene,
    :Sucrose,
    :Terphenyl,
    :Testes_Icrp,
    :Tetrachloroethylene,
    :Thallium_Chloride,
    :Tissue_Soft_Icrp,
    :Tissue_Soft_Icru_Four_Component,
    :Tissue_Equivalent_Gas_Methane_Based,
    :Tissue_Equivalent_Gas_Propane_Based,
    :Titanium_Dioxide,
    :Toluene,
    :Trichloroethylene,
    :Triethyl_Phosphate,
    :Tungsten_Hexafluoride,
    :Uranium_Dicarbide,
    :Uranium_Monocarbide,
    :Uranium_Oxide,
    :Urea,
    :Valine,
    :Viton_Fluoroelastomer,
    :Water,
    :Water_Vapor,
    :Xylene,
    :Carbon_Graphite,
]
