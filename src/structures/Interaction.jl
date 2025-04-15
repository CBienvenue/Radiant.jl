
abstract type Interaction end

# Method(s)
"""
    get_in_particles(interaction::Interaction)

Get the interaction incoming particle.

# Input Argument(s)
- `interaction::Interaction` : interaction.

# Output Argument(s)
- `incoming_particle::Type` : incoming particle.

"""
function get_in_particles(interaction::Interaction)
    return interaction.incoming_particle
end

"""
    get_out_particles(interaction::Interaction)

Get the interaction outgoing particles.

# Input Argument(s)
- `interaction::Interaction` : interaction.

# Output Argument(s)
- `interaction_particles::Vector{Type}` : outgoing particle.

"""
function get_out_particles(interaction::Interaction)
    return interaction.interaction_particles
end

"""
    get_is_subshells_dependant(interaction::Interaction)

Is the interaction cross-sections are subshells dependant.

# Input Argument(s)
- `interaction::Interaction` : interaction.

# Output Argument(s)
- `is_subshells_dependant::Bool` : is the interaction cross-sections are subshells
  dependant.

"""
function get_is_subshells_dependant(interaction::Interaction)
    return interaction.is_subshells_dependant
end


"""
    get_types(interaction::Interaction,particle_in::Type,particle_out::Type)

Get the interaction types.

# Input Argument(s)
- `interaction::Interaction` : interaction.
- `particle_in::Type` : incoming particle.
- `particle_out::Type` : outgoing particle.

# Output Argument(s)
- `types::Vector{String}` : types of interaction.

"""
function get_types(interaction::Interaction,particle_in::Type,particle_out::Type)
    if haskey(interaction.interaction_types,(particle_in,particle_out))
        return interaction.interaction_types[(particle_in,particle_out)]
    else
        return Vector{String}()
    end
end

"""
    get_mass_energy_variation(interaction::Interaction,type::String="")

Get the mass energy variation in the interaction.

# Input Argument(s)
- `interaction::Interaction` : interaction.
- `type::String` : type of interaction.
- `is_production::Bool` : is the mass energy variation comes from production interaction.

# Output Argument(s)
- `ΔQ::Float64` : types of interaction.

"""
function get_mass_energy_variation(interaction::Interaction,type::String,is_production::Bool)
    if typeof(interaction) == Pair_Production && ~is_production
        return 2.0 # One electron and one positron are created
    elseif typeof(interaction) == Annihilation
        if type ∈ ["P_inel","P_brems"] && is_production
            return -1.0 # One electron and one positron are annihilated
        elseif type ∈ ["P_pp"] && is_production
            return 0.0 # One electron and one positron are created, then the positron annihilate with an electron
        elseif ~is_production
            return -2.0
        else
            return 0.0
        end
    else
        return 0.0 # No particle are created nor annihilated
    end
end

"""
    get_charge_variation(interaction::Interaction,type::String)

Get the charge deposited and extracted in the interaction.

# Input Argument(s)
- `interaction::Interaction` : interaction.
- `type::String` : type of interaction.
- `charge_in::Real` : type of interaction.
- `charge_out::Real` : type of interaction.

# Output Argument(s)
- `q_deposited::Float64` : charge deposited by the incoming particle.
- `q_extracted::Float64` : charge extracted by the production of particle(s).

"""
function get_charge_variation(interaction::Interaction,type::String,charge_in::Real,charge_out::Real)
    if typeof(interaction) == Annihilation
        if type == "P_pp"
            q_deposited = 0; q_extracted = 1/2
        else
            q_deposited = 1; q_extracted = 0
        end
    elseif typeof(interaction) == Pair_Production
        q_deposited = 0; q_extracted = 0 # no electron extracted from medium in pair production
    else
        q_deposited = -charge_in; q_extracted = -charge_out
    end
    return q_deposited, q_extracted
end

"""
    get_scattering_model(interaction::Interaction)

Get the scattering model associated with the interaction.

# Input Argument(s)
- `interaction::Interaction` : interaction.

# Output Argument(s)
- `scattering_model::String` : scattering model.
- `is_CSD::Bool` : boolean indicating if continuous slowing-down term are employed.
- `is_AFP::Bool` : boolean indicating if angular Fokker-Planck term are employed.
- `is_AFP_decomposition::Bool` : boolean indicating if angular Fokker-Planck decomposition
  method is employed.

"""
function get_scattering_model(interaction::Interaction)

    if hasfield(typeof(interaction), :is_ETC)
        is_ETC = interaction.is_ETC
    else
        is_ETC = false
    end
    
    return interaction.scattering_model, interaction.is_CSD, interaction.is_AFP, interaction.is_AFP_decomposition, is_ETC
end

"""
    is_elastic_scattering(interaction::Interaction)

Indicate if the scattering is elastic or not.

# Input Argument(s)
- `interaction::Interaction` : interaction.

# Output Argument(s)
- `is_elastic::Bool` : boolean indicating if the scattering is elastic.

"""
function is_elastic_scattering(interaction::Interaction)
    return interaction.is_elastic
end

"""
    in_distribution_dispatch(interaction::Interaction)

Describe the energy discretization method for the incoming particle.

# Input Argument(s)
- `this::Interaction` : Interaction structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution_dispatch(interaction::Interaction)
    itype = typeof(interaction)
    if itype == Annihilation
        return in_distribution(interaction)
    elseif itype == Bremsstrahlung
        return in_distribution(interaction)
    elseif itype == Compton
        return in_distribution(interaction)
    elseif itype == Elastic_Collision
        return in_distribution(interaction)
    elseif itype == Inelastic_Collision
        return in_distribution(interaction)
    elseif itype == Pair_Production
        return in_distribution(interaction)
    elseif itype == Photoelectric
        return in_distribution(interaction)
    elseif itype == Rayleigh
        return in_distribution(interaction)
    elseif itype == Relaxation
        return in_distribution(interaction)
    else
        error("Unknown interaction.")
    end
end

"""
    out_distribution_dispatch(interaction::Interaction,type::String)

Describe the energy discretization method for the outgoing particle.

# Input Argument(s)
- `this::Interaction` : Interaction structure.
- `type::String` : type of interaction.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution_dispatch(interaction::Interaction,type::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return out_distribution(interaction,type)
    elseif itype == Bremsstrahlung
        return out_distribution(interaction)
    elseif itype == Compton
        return out_distribution(interaction)
    elseif itype == Elastic_Collision
        return out_distribution(interaction)
    elseif itype == Inelastic_Collision
        return out_distribution(interaction)
    elseif itype == Pair_Production
        return out_distribution(interaction)
    elseif itype == Photoelectric
        return out_distribution(interaction,type)
    elseif itype == Rayleigh
        return out_distribution(interaction)
    elseif itype == Relaxation
        return out_distribution(interaction)
    else
        error("Unknown interaction.")
    end
end

"""
    bounds_dispatch(interaction::Interaction,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,
    gi::Int64,gf::Int64,type::String,Ui::Float64,Ec::Float64,particle::Particle)

Gives the integration energy bounds for the outgoing particle. 

# Input Argument(s)
- `this::Interaction` : Interaction structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `gi::Int64` : group of the incoming particle.
- `gf::Int64` : group of the outgoing particle.
- `type::String` : type of interaction.
- `Ui::Float64` : binding energy of the electron in the subshell.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds_dispatch(interaction::Interaction,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,gi::Int64,gf::Int64,type::String,Ui::Float64,Ec::Float64,particle::Particle)
    itype = typeof(interaction)
    if itype == Annihilation
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Bremsstrahlung
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type,Ec)
    elseif itype == Compton
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Elastic_Collision
        return bounds(interaction,Ef⁻,Ef⁺,gi,gf)
    elseif itype == Inelastic_Collision
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type,Ec,Ui,particle)
    elseif itype == Pair_Production
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Photoelectric
        return bounds(interaction,Ef⁻,Ef⁺)
    elseif itype == Rayleigh
        return bounds(interaction,Ef⁻,Ef⁺,gi,gf)
    elseif itype == Relaxation
        return bounds(interaction,Ef⁻,Ef⁺)
    else
        error("Unknown interaction.")
    end
end

"""
    dcs_dispatch(interaction::Interaction,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,
    particle::Particle,type::String,iz::Int64,particles::Vector{Particle},
    Ein::Vector{Float64},vec_Z::Vector{Int64},Ef⁻::Float64,Ef⁺::Float64,δi::Int64,
    Ui::Float64,Zi::Real,Ti::Float64,Ec::Float64,incoming_particle::Particle)

Gives the Legendre moments of the scattering cross-sections. 

# Input Argument(s)
- `this::Interaction` : Interaction structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `Z::Int64` : atomic number.
- `particle::Particle` : scattered or produced particle.
- `type::String` : type of interaction.
- `iz::Int64` : index of the element in the material.
- `particles::Vector{Particle}` : list of particles.
- `Ein::Vector{Float64}` : energy boundaries associated with the incoming particle.
- `vec_Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Ef⁻::Float64` : upper bound of the outgoing particle group.
- `Ef⁺::Float64` : lower bound of the outgoing particle group.
- `δi::Int64` : subshell index.
- `Ui::Float64` : binding energy of the electron on the subshell.
- `Zi::Real` : mean number of electron on the subshell.
- `Ti::Float64` : mean kinetic energy of the electron on the subshell.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `incoming_particle::Particle` : incoming particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs_dispatch(interaction::Interaction,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,scattered_particle::Particle,type::String,iz::Int64,particles::Vector{Particle},Ein::Vector{Float64},vec_Z::Vector{Int64},Ef⁻::Float64,Ef⁺::Float64,δi::Int64,Ui::Float64,Zi::Real,Ti::Float64,Ec::Float64,incoming_particle::Particle,ρ::Float64)
    itype = typeof(interaction)
    if itype == Annihilation
        return dcs(interaction,L,Ei,Ef,type,Z,Ein,Ec)
    elseif itype == Bremsstrahlung
        return dcs(interaction,L,Ei,Ef,Z,incoming_particle,type,iz)
    elseif itype == Compton
        return dcs(interaction,L,Ei,Ef,type,Z,iz,δi,scattered_particle)
    elseif itype == Elastic_Collision
        return dcs(interaction,L,Ei,Z,incoming_particle,Ein[end])
    elseif itype == Inelastic_Collision
        return dcs(interaction,L,Ei,Ef,type,incoming_particle,Ui,Zi,Ti)
    elseif itype == Pair_Production
        return dcs(interaction,L,Ei,Ef,Z,type,iz,particles)
    elseif itype == Photoelectric
        return dcs(interaction,L,Ei,Z,δi,Ef⁻,Ef⁺)
    elseif itype == Rayleigh
        return dcs(interaction,L,Ei,Z,iz)
    elseif itype == Relaxation
        return dcs(interaction,L,Ei,Ein[end],Ec,Z,δi,Ef⁻,Ef⁺,incoming_particle,scattered_particle)
    else
        error("Unknown interaction.")
    end
end

"""
    tcs_dispatch(interaction::Interaction,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,
    particle::Particle,Ecutoff::Float64,Eout::Vector{Float64},type::String)

Gives the total cross-section. 

# Input Argument(s)
- `this::Interaction` : Interaction structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `iz::Int64` : index of the element in the material.
- `particle::Particle` : incoming particle.
- `Ecutoff::Float64` : cutoff energy.
- `Eout::Vector{Float64}` : outgoing energy boundaries.
- `type::String` : type of interaction.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs_dispatch(interaction::Interaction,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::Particle,Ecutoff::Float64,Eout::Vector{Float64},ρ::Float64)
    itype = typeof(interaction)
    if itype == Annihilation
        return tcs(interaction,Ei,Z)
    elseif itype == Bremsstrahlung
        return tcs(interaction,Ei,Z,Ec,particle,Eout)
    elseif itype == Compton
        return tcs(interaction,Ei,Z,Eout)
    elseif itype == Elastic_Collision
        return tcs(interaction,Ei,Z,particle,Ecutoff)
    elseif itype == Inelastic_Collision
        return tcs(interaction,Ei,Ec,particle,Z)
    elseif itype == Pair_Production
        return tcs(interaction,Ei,Z,Eout)
    elseif itype == Photoelectric
        return tcs(interaction,Ei,Z)
    elseif itype == Rayleigh
        return tcs(interaction,Ei,Z,iz)
    else
        error("Unknown interaction.")
    end
end

"""
    acs_dispatch(interaction::Interaction,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,
    particle::Particle,Ecutoff::Float64,Eout::Vector{Float64},type::String)

Gives the absorption cross-section. 

# Input Argument(s)
- `this::Interaction` : Interaction structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `iz::Int64` : index of the element in the material.
- `particle::Particle` : incoming particle.
- `Ecutoff::Float64` : cutoff energy.
- `Eout::Vector{Float64}` : outgoing energy boundaries.
- `type::String` : type of interaction.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs_dispatch(interaction::Interaction,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::Particle,Ecutoff::Float64,Eout::Vector{Float64},ρ::Float64)
    itype = typeof(interaction)
    if itype == Annihilation
        return acs(interaction,Ei,Z)
    elseif itype == Bremsstrahlung
        return acs(interaction,Ei,Z,Ec,particle,Ecutoff)
    elseif itype == Compton
        return acs(interaction,Ei,Z,Ecutoff)
    elseif itype == Elastic_Collision
        return acs(interaction)
    elseif itype == Inelastic_Collision
        return acs(interaction,Ei,Ec,particle,Z,Ecutoff)
    elseif itype == Pair_Production
        return acs(interaction,Ei,Z,Eout)
    elseif itype == Photoelectric
        return acs(interaction,Ei,Z)
    elseif itype == Rayleigh
        return acs(interaction)
    else
        error("Unknown interaction.")
    end
end

"""
    sp_dispatch(interaction::Interaction,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,
    state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle,
    Eout::Vector{Float64})

Gives the stopping power.

# Input Argument(s)
- `this::Interaction` : Interaction structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `state_of_matter::String` : state of matter.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp_dispatch(interaction::Interaction,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle,Eout::Vector{Float64})
    itype = typeof(interaction)
    if itype == Bremsstrahlung
        return sp(interaction,Z,ωz,ρ,Ei,Ec,Eout,particle)
    elseif itype == Inelastic_Collision
        return sp(interaction,Z,ωz,ρ,state_of_matter,Ei,Ec,particle)
    else
        error("Unknown interaction.")
    end
end

"""
    mt_dispatch(interaction::Interaction)

Gives the momentum transfer.

# Input Argument(s)
- `this::Interaction` : Interaction structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt_dispatch(interaction::Interaction)
    itype = typeof(interaction)
    if itype == Bremsstrahlung
        return mt(interaction)
    elseif itype == Inelastic_Collision
        return mt(interaction)
    else
        error("Unknown interaction.")
    end
end