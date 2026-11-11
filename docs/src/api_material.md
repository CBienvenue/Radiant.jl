# Material

## Structure
```@docs
Radiant.Material
```

## Constructors and Setters
```@docs
Radiant.set_tag(this::Radiant.Material,tag::String)
Radiant.set_density(this::Radiant.Material,density::Real)
Radiant.set_state_of_matter(this::Radiant.Material,state_of_matter::String)
Radiant.add_element(this::Radiant.Material,symbol::String,weight_fraction::Real)
Radiant.set_mean_excitation_energy(this::Radiant.Material,I::Real)
Radiant.set_stopping_power(this::Radiant.Material,particle::Type,file::String)
```

## Getters
```@docs
Radiant.get_tag(this::Radiant.Material)
Radiant.get_density(this::Radiant.Material)
Radiant.get_state_of_matter(this::Radiant.Material)
Radiant.get_number_of_elements(this::Radiant.Material)
Radiant.get_atomic_numbers(this::Radiant.Material)
Radiant.get_weight_fractions(this::Radiant.Material)
Radiant.get_mean_excitation_energy(this::Radiant.Material)
Radiant.get_stopping_power(this::Radiant.Material,particle::Type)
```

## Printing
```@docs
Radiant.println(this::Radiant.Material)
```
