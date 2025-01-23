# 2 Materials

Materials in Radiant are defined by a `Radiant.Material` object. Once a material has been instantiated, the elemental composition of the material can be defined with its density. For example, the following code would be use to define water:

```julia
water = Material()
water.set_density(1.0)           # Density of 1 g/cm³
water.add_element("H",0.1111)    # Hydrogen with mass weight of 11.11 %
water.add_element("O",0.8889)    # Oxygen with mass weight of 88.89 %
```

and for a monoelemental element, like aluminium, the following can be used:

```julia
al = Material()
al.add_element("Al")             # Aluminium with mass weight of 100 %
```

Note that the absence of explicit mass weight means the material is monoelemental and the absence of explicit density correspond to the material density at 20°C. The material information can be printed using `println()` function. For example, for the object `water` created before, the function `println(water)` generates the following output:
```
Material:
   ID:                           53
   Density (g/cm³):              1.0
   Elements in the compound:     ["h", "o"]
   Weight fractions:             [0.1111, 0.8889]
   State of matter:              liquid
```
where the ID field correspond to an internal ID to uniquely identify the object. An empty field can take the value `missing` if the data is required for calculations, and simply `nothing` if not.