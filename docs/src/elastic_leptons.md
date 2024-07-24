## Structure
```@docs
Radiant.Elastic_Leptons
```

## Methods
```@docs
Radiant.set_model(this::Radiant.Elastic_Leptons,model::String)
Radiant.set_transport_correction(this::Radiant.Elastic_Leptons,is_ETC::Bool)
Radiant.set_angular_fokker_planck(this::Radiant.Elastic_Leptons,is_AFP::Bool)
Radiant.set_kawrakow_correction(this::Radiant.Elastic_Leptons,is_kawrakow_correction::Bool,subshell_dependant_inelastic::Bool=true)
Radiant.set_interaction_types(this::Radiant.Elastic_Leptons,interaction_types::Dict{Tuple{String,String},Vector{String}})
```