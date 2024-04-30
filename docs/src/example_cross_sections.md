# Production of cross-sections

RADIANT can produce multigroup cross-sections for photons, electrons and positrons.

## Example 1: Production of coupled photons-electrons-positrons cross-sections in water

```julia

# Define the medium (water)
water = Material()
water.set_name("water")
water.set_density(1.0)
water.set_state_of_matter("liquid")
water.add_element("H",0.1111)
water.add_element("O",0.8889)

# Define the cross-sections parameters
cs = Cross_Sections()
cs.set_source("RADIANT")
cs.set_materials([water]) 
cs.set_particles(["photons","electrons","positrons"])
cs.set_energy(10)
cs.set_cutoff(0.001)
cs.set_number_of_groups([80,80,80])
cs.set_group_structure(["log","log","log"])
cs.set_solvers(["BTE","BFP","BFP"])
interaction_list = [
    Elastic_Leptons(),
    Inelastic_Leptons(),
    Bremsstrahlung(),
    Annihilation(),
    Auger(),
    Fluorescence(),
    Pair_Production(),
    Photoelectric(),
    Rayleigh(),
    Compton()
]
cs.set_interactions(interaction_list)
cs.set_legendre_order(7)
cs.build()

# Save cross-sections in a FMAC-M file
cs.write("(path here)/fmac_m.txt")

```

