function interaction_interdependances(interactions::Vector{Interaction},particles::Vector{Particle})

# Loop over each interaction
for interaction_i in interactions 

    #----
    # Annihilation
    #----
    if typeof(interaction_i) == Annihilation && any(is_photon.(particles))

        # Annihilation of positrons scattered under the cutoff from inelastic collisionnal interaction
        if any(is_positron.(particles)) && haskey(interaction_i.interaction_types,(Positron,Photon)) && "P_inel" ∈ interaction_i.interaction_types[(Positron,Photon)]
            is_interaction = false
            for interaction_j in interactions 
                if typeof(interaction_j) == Inelastic_Collision
                    interaction_i.inelastic_collision_model = interaction_j
                    is_interaction = true
                    break
                end
            end
            if ~is_interaction
                filter!(x -> x != "P_inel", interaction_i.interaction_types[(Positron,Photon)])
            end
        end

        # Annihilation of positrons scattered under the cutoff from Bremsstrahlung interaction
        if any(is_positron.(particles)) && haskey(interaction_i.interaction_types,(Positron,Photon)) && "P_brems" ∈ interaction_i.interaction_types[(Positron,Photon)]
            is_interaction = false
            for interaction_j in interactions 
                if typeof(interaction_j) == Bremsstrahlung
                    interaction_i.bremsstrahlung_model = interaction_j
                    is_interaction = true
                    break
                end
            end
            if ~is_interaction
                filter!(x -> x != "P_brems", interaction_i.interaction_types[(Positron,Photon)])
            end
        end

        # Annihilation of positrons produced under the cutoff following pair production event
        if haskey(interaction_i.interaction_types,(Photon,Photon)) && "P_pp" ∈ interaction_i.interaction_types[(Photon,Photon)]
            is_interaction = false
            for interaction_j in interactions 
                if typeof(interaction_j) == Pair_Production
                    interaction_i.pair_production_model = interaction_j
                    is_interaction = true
                    break
                end
            end
            if ~is_interaction
                filter!(x -> x != "P_pp", interaction_i.interaction_types[(Photon,Photon)])
            end
        end
    end

    #----
    # Elastic Collision
    #----
    if typeof(interaction_i) == Elastic_Collision

        # Model to use for Kawrakow correction
        if interaction_i.is_kawrakow_correction
            for interaction_j in interactions
                if typeof(interaction_j) == Inelastic_Collision
                    interaction_i.is_subshell_inelastic = interaction_j.is_subshells_dependant
                    break
                end
            end
        end
    end

    #----
    # Pair production
    #----
    if typeof(interaction_i) == Pair_Production

        # Generate two photons, if there is not explicit positron transport
        if ~any(is_positron.(particles))
            interaction_i.is_positron = false
        end
    end
    
end
end