#----------------------------------------------------------------#
#
# This file defines functions to create materials corresponding to the NIST database of elements. Each function creates a `Material` object, sets its density and mean excitation energy, and adds the appropriate element with its weight fraction.
#
# References:
# - National Institute of Standards and Technology (NIST) database: https://physics.nist.gov/PhysRefData/Star/Text/download.html.
#
#----------------------------------------------------------------#

"""
    Hydrogen()

NIST material: HYDROGEN
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

NIST material: HELIUM
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

NIST material: LITHIUM
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

NIST material: BERYLLIUM
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

NIST material: BORON
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

NIST material: AMORPHOUS CARBON
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

NIST material: NITROGEN
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

NIST material: OXYGEN
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

NIST material: FLUORINE
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

NIST material: NEON
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

NIST material: SODIUM
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

NIST material: MAGNESIUM
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

NIST material: ALUMINUM
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

NIST material: SILICON
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

NIST material: PHOSPHORUS
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

NIST material: SULFUR
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

NIST material: CHLORINE
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

NIST material: ARGON
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

NIST material: POTASSIUM
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

NIST material: CALCIUM
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

NIST material: SCANDIUM
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

NIST material: TITANIUM
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

NIST material: VANADIUM
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

NIST material: CHROMIUM
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

NIST material: MANGANESE
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

NIST material: IRON
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

NIST material: COBALT
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

NIST material: NICKEL
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

NIST material: COPPER
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

NIST material: ZINC
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

NIST material: GALLIUM
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

NIST material: GERMANIUM
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

NIST material: ARSENIC
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

NIST material: SELENIUM
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

NIST material: BROMINE
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

NIST material: KRYPTON
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

NIST material: RUBIDIUM
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

NIST material: STRONTIUM
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

NIST material: YTTRIUM
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

NIST material: ZIRCONIUM
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

NIST material: NIOBIUM
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

NIST material: MOLYBDENUM
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

NIST material: TECHNETIUM
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

NIST material: RUTHENIUM
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

NIST material: RHODIUM
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

NIST material: PALLADIUM
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

NIST material: SILVER
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

NIST material: CADMIUM
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

NIST material: INDIUM
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

NIST material: TIN
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

NIST material: ANTIMONY
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

NIST material: TELLURIUM
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

NIST material: IODINE
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

NIST material: XENON
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

NIST material: CESIUM
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

NIST material: BARIUM
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

NIST material: LANTHANUM
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

NIST material: CERIUM
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

NIST material: PRASEODYMIUM
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

NIST material: NEODYMIUM
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

NIST material: PROMETHIUM
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

NIST material: SAMARIUM
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

NIST material: EUROPIUM
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

NIST material: GADOLINIUM
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

NIST material: TERBIUM
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

NIST material: DYSPROSIUM
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

NIST material: HOLMIUM
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

NIST material: ERBIUM
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

NIST material: THULIUM
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

NIST material: YTTERBIUM
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

NIST material: LUTETIUM
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

NIST material: HAFNIUM
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

NIST material: TANTALUM
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

NIST material: TUNGSTEN
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

NIST material: RHENIUM
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

NIST material: OSMIUM
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

NIST material: IRIDIUM
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

NIST material: PLATINUM
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

NIST material: GOLD
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

NIST material: MERCURY
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

NIST material: THALLIUM
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

NIST material: LEAD
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

NIST material: BISMUTH
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

NIST material: POLONIUM
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

NIST material: ASTATINE
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

NIST material: RADON
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

NIST material: FRANCIUM
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

NIST material: RADIUM
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

NIST material: ACTINIUM
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

NIST material: THORIUM
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

NIST material: PROTACTINIUM
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

NIST material: URANIUM
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

NIST material: NEPTUNIUM
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

NIST material: PLUTONIUM
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

NIST material: AMERICIUM
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

NIST material: CURIUM
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

NIST material: BERKELIUM
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

NIST material: CALIFORNIUM
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

NIST material: A-150 TISSUE-EQUIVALENT PLASTIC
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

NIST material: ACETONE
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

NIST material: ACETYLENE
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

NIST material: ADENINE
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

NIST material: ADIPOSE TISSUE (ICRP)
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

NIST material: AIR, DRY (NEAR SEA LEVEL)
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

NIST material: ALANINE
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

NIST material: ALUMINUM OXIDE
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

NIST material: AMBER
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

NIST material: AMMONIA
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

NIST material: ANILINE
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

NIST material: ANTHRACENE
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

NIST material: B-100 BONE-EQUIVALENT PLASTIC
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

NIST material: BAKELITE
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

NIST material: BARIUM FLUORIDE
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

NIST material: BARIUM SULFATE
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

NIST material: BENZENE
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

NIST material: BERYLLIUM OXIDE
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

NIST material: BISMUTH GERMANIUM OXIDE
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

NIST material: BLOOD (ICRP)
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

NIST material: BONE, COMPACT (ICRU)
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

NIST material: BONE, CORTICAL (ICRP)
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

NIST material: BORON CARBIDE
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

NIST material: BORON OXIDE
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

NIST material: BRAIN (ICRP)
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

NIST material: BUTANE
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

NIST material: N-BUTYL ALCOHOL
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

NIST material: C-552 AIR-EQUIVALENT PLASTIC
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

NIST material: CADMIUM TELLURIDE
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

NIST material: CADMIUM TUNGSTATE
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

NIST material: CALCIUM CARBONATE
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

NIST material: CALCIUM FLUORIDE
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

NIST material: CALCIUM OXIDE
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

NIST material: CALCIUM SULFATE
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

NIST material: CALCIUM TUNGSTATE
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

NIST material: CARBON DIOXIDE
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

NIST material: CARBON TETRACHLORIDE
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

NIST material: CELLULOSE ACETATE, CELLOPHANE
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

NIST material: CELLULOSE ACETATE BUTYRATE
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

NIST material: CELLULOSE NITRATE
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

NIST material: CERIC SULFATE DOSIMETER SOLUTION
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

NIST material: CESIUM FLUORIDE
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

NIST material: CESIUM IODIDE
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

NIST material: CHLOROBENZENE
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

NIST material: CHLOROFORM
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

NIST material: CONCRETE, PORTLAND
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

NIST material: CYCLOHEXANE
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

NIST material: 1,2-DICHLOROBENZENE
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

NIST material: DICHLORODIETHYL ETHER
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

NIST material: 1,2-DICHLOROETHANE
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

NIST material: DIETHYL ETHER
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

NIST material: N,N-DIMETHYL FORMAMIDE
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

NIST material: DIMETHYL SULFOXIDE
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

NIST material: ETHANE
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

NIST material: ETHYL ALCOHOL
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

NIST material: ETHYL CELLULOSE
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

NIST material: ETHYLENE
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

NIST material: EYE LENS (ICRP)
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

NIST material: FERRIC OXIDE
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

NIST material: FERROBORIDE
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

NIST material: FERROUS OXIDE
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

NIST material: FERROUS SULFATE DOSIMETER SOLUTION
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

NIST material: FREON-12
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

NIST material: FREON-12B2
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

NIST material: FREON-13
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

NIST material: FREON-13B1
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

NIST material: FREON-13I1
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

NIST material: GADOLINIUM OXYSULFIDE
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

NIST material: GALLIUM ARSENIDE
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

NIST material: GEL IN PHOTOGRAPHIC EMULSION
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

NIST material: Pyrex Glass
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

NIST material: GLASS, LEAD
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

NIST material: GLASS, PLATE
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

NIST material: GLUCOSE
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

NIST material: GLUTAMINE
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

NIST material: GLYCEROL
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

NIST material: GUANINE
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

NIST material: GYPSUM, PLASTER OF PARIS
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

NIST material: N-HEPTANE
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

NIST material: N-HEXANE
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

NIST material: KAPTON POLYIMIDE FILM
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

NIST material: LANTHANUM OXYBROMIDE
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

NIST material: LANTHANUM OXYSULFIDE
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

NIST material: LEAD OXIDE
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

NIST material: LITHIUM AMIDE
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

NIST material: LITHIUM CARBONATE
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

NIST material: LITHIUM FLUORIDE
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

NIST material: LITHIUM HYDRIDE
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

NIST material: LITHIUM IODIDE
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

NIST material: LITHIUM OXIDE
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

NIST material: LITHIUM TETRABORATE
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

NIST material: LUNG (ICRP)
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

NIST material: M3 WAX
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

NIST material: MAGNESIUM CARBONATE
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

NIST material: MAGNESIUM FLUORIDE
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

NIST material: MAGNESIUM OXIDE
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

NIST material: MAGNESIUM TETRABORATE
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

NIST material: MERCURIC IODIDE
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

NIST material: METHANE
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

NIST material: METHANOL
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

NIST material: MIX D WAX
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

NIST material: MS20 TISSUE SUBSTITUTE
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

NIST material: MUSCLE, SKELETAL (ICRP)
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

NIST material: MUSCLE, STRIATED (ICRU)
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

NIST material: MUSCLE-EQUIVALENT LIQUID, WITH SUCROSE
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

NIST material: MUSCLE-EQUIVALENT LIQUID, WITHOUT SUCROSE
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

NIST material: NAPHTHALENE
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

NIST material: NITROBENZENE
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

NIST material: NITROUS OXIDE
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

NIST material: NYLON, DU PONT ELVAMIDE 8062
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

NIST material: NYLON, TYPE 6 AND TYPE 6/6
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

NIST material: NYLON, TYPE 6/10
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

NIST material: NYLON, TYPE 11 (RILSAN)
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

NIST material: OCTANE, LIQUID
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

NIST material: PARAFFIN WAX
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

NIST material: N-PENTANE
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

NIST material: PHOTOGRAPHIC EMULSION
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

NIST material: PLASTIC SCINTILLATOR (VINYLTOLUENE BASED)
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

NIST material: PLUTONIUM DIOXIDE
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

NIST material: POLYACRYLONITRILE
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

NIST material: POLYCARBONATE (MAKROLON, LEXAN)
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

NIST material: POLYCHLOROSTYRENE
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

NIST material: POLYETHYLENE
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

NIST material: POLYETHYLENE TEREPHTHALATE (MYLAR)
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

NIST material: POLYMETHYL METHACRALATE (LUCITE, PERSPEX, PLEXIGLASS)
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

NIST material: POLYOXYMETHYLENE
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

NIST material: POLYPROPYLENE
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

NIST material: POLYSTYRENE
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

NIST material: POLYTETRAFLUOROETHYLENE (TEFLON)
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

NIST material: POLYTRIFLUOROCHLOROETHYLENE
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

NIST material: POLYVINYL ACETATE
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

NIST material: POLYVINYL ALCOHOL
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

NIST material: POLYVINYL BUTYRAL
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

NIST material: POLYVINYL CHLORIDE
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

NIST material: POLYVINYLIDENE CHLORIDE, SARAN
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

NIST material: POLYVINYLIDENE FLUORIDE
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

NIST material: POLYVINYL PYRROLIDONE
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

NIST material: POTASSIUM IODIDE
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

NIST material: POTASSIUM OXIDE
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

NIST material: PROPANE
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

NIST material: PROPANE, LIQUID
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

NIST material: N-PROPYL ALCOHOL
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

NIST material: PYRIDINE
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

NIST material: RUBBER, BUTYL
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

NIST material: RUBBER, NATURAL
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

NIST material: RUBBER, NEOPRENE
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

NIST material: SILICON DIOXIDE
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

NIST material: SILVER BROMIDE
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

NIST material: SILVER CHLORIDE
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

NIST material: SILVER HALIDES IN PHOTOGRAPHIC EMULSION
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

NIST material: SILVER IODIDE
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

NIST material: SKIN (ICRP)
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

NIST material: SODIUM CARBONATE
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

NIST material: SODIUM IODIDE
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

NIST material: SODIUM MONOXIDE
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

NIST material: SODIUM NITRATE
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

NIST material: STILBENE
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

NIST material: SUCROSE
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

NIST material: TERPHENYL
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

NIST material: TESTES (ICRP)
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

NIST material: TETRACHLOROETHYLENE
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

NIST material: THALLIUM CHLORIDE
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

NIST material: TISSUE, SOFT (ICRP)
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

NIST material: TISSUE, SOFT (ICRU FOUR-COMPONENT)
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

NIST material: TISSUE-EQUIVALENT GAS (METHANE BASED)
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

NIST material: TISSUE-EQUIVALENT GAS (PROPANE BASED)
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

NIST material: TITANIUM DIOXIDE
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

NIST material: TOLUENE
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

NIST material: TRICHLOROETHYLENE
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

NIST material: TRIETHYL PHOSPHATE
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

NIST material: TUNGSTEN HEXAFLUORIDE
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

NIST material: URANIUM DICARBIDE
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

NIST material: URANIUM MONOCARBIDE
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

NIST material: URANIUM OXIDE
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

NIST material: UREA
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

NIST material: VALINE
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

NIST material: VITON FLUOROELASTOMER
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

NIST material: WATER, LIQUID
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

NIST material: WATER VAPOR
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

NIST material: XYLENE
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

NIST material: CARBON (GRAPHITE)
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
