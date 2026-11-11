# Inelastic_Collision

## Structure
```@docs
Radiant.Inelastic_Collision
```

## Setters
```@docs
Radiant.set_interaction_types(this::Radiant.Inelastic_Collision,interaction_types)
Radiant.set_density_correction(this::Radiant.Inelastic_Collision,density_correction::String)
Radiant.set_is_shell_correction(this::Radiant.Inelastic_Collision,is_shell_correction::Bool)
Radiant.set_is_subshells_dependant(this::Radiant.Inelastic_Collision,is_subshells_dependant::Bool)
Radiant.set_scattering_model(this::Radiant.Inelastic_Collision,scattering_model::String)
Radiant.set_ecut(this::Radiant.Inelastic_Collision,particle::Type,energy::Real)
Radiant.set_emax(this::Radiant.Inelastic_Collision,particle::Type,energy::Real)
Radiant.set_low_energy_model(this::Radiant.Inelastic_Collision,model::String)
Radiant.set_low_energy_source(this::Radiant.Inelastic_Collision,source::String)
Radiant.set_is_bragg_additivity(this::Radiant.Inelastic_Collision,is_bragg_additivity::Bool)
Radiant.set_use_srim_sp(this::Radiant.Inelastic_Collision,use_srim_sp::Bool)
```

## Getters
```@docs
Radiant.get_ecut(this::Radiant.Inelastic_Collision,particle::Type)
Radiant.get_emax(this::Radiant.Inelastic_Collision,particle::Type)
Radiant.get_low_energy_model(this::Radiant.Inelastic_Collision)
Radiant.get_low_energy_source(this::Radiant.Inelastic_Collision)
Radiant.get_is_bragg_additivity(this::Radiant.Inelastic_Collision)
Radiant.get_use_srim_sp(this::Radiant.Inelastic_Collision)
```
