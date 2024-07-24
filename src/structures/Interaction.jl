
abstract type Interaction end

# Method(s)
function get_in_particles(interaction::Interaction)
    return interaction.incoming_particle
end

function get_out_particles(interaction::Interaction)
    return interaction.interaction_particles
end

function get_types(interaction::Interaction,particle_in::String,particle_out::String)
    if haskey(interaction.interaction_types,(particle_in,particle_out))
        return interaction.interaction_types[(particle_in,particle_out)]
    else
        return Vector{String}()
    end
end

function in_distribution_dispatch(interaction::Interaction,type::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return in_distribution(interaction)
    elseif itype == Bremsstrahlung
        return in_distribution(interaction)
    elseif itype == Compton
        return in_distribution(interaction)
    elseif itype == Elastic_Leptons
        return in_distribution(interaction)
    elseif itype == Inelastic_Leptons
        return in_distribution(interaction)
    elseif itype == Pair_Production
        return in_distribution(interaction)
    elseif itype == Photoelectric
        return in_distribution(interaction)
    elseif itype == Rayleigh
        return in_distribution(interaction)
    elseif itype == Fluorescence
        return in_distribution(interaction)
    elseif itype == Auger
        return in_distribution(interaction)
    else
        error("Unknown interaction.")
    end
end

function out_distribution_dispatch(interaction::Interaction,type::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return out_distribution(interaction,type)
    elseif itype == Bremsstrahlung
        return out_distribution(interaction)
    elseif itype == Compton
        return out_distribution(interaction)
    elseif itype == Elastic_Leptons
        return out_distribution(interaction)
    elseif itype == Inelastic_Leptons
        return out_distribution(interaction)
    elseif itype == Pair_Production
        return out_distribution(interaction)
    elseif itype == Photoelectric
        return out_distribution(interaction,type)
    elseif itype == Rayleigh
        return out_distribution(interaction)
    elseif itype == Fluorescence
        return out_distribution(interaction)
    elseif itype == Auger
        return out_distribution(interaction)
    else
        error("Unknown interaction.")
    end
end

function bounds_dispatch(interaction::Interaction,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,gi::Int64,gf::Int64,type::String,Ui::Float64,Z::Int64,Ein::Vector{Float64},Ngi::Int64,Ec::Float64,particle::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Bremsstrahlung
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type,Ec)
    elseif itype == Compton
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Elastic_Leptons
        return bounds(interaction,Ef⁻,Ef⁺,gi,gf)
    elseif itype == Inelastic_Leptons
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type,Ec,Ui,particle)
    elseif itype == Pair_Production
        return bounds(interaction,Ef⁻,Ef⁺,Ei,type)
    elseif itype == Photoelectric
        return bounds(interaction,Ef⁻,Ef⁺,gi,type,Ui,Ein)
    elseif itype == Rayleigh
        return bounds(interaction,Ef⁻,Ef⁺,gi,gf)
    elseif itype == Fluorescence
        return bounds(interaction,Ef⁻,Ef⁺,gi,type,Ui,Ein)
    elseif itype == Auger
        return bounds(interaction,Ef⁻,Ef⁺,gi,type,Ui,Ein)
    else
        error("Unknown interaction.")
    end
end

function dcs_dispatch(interaction::Interaction,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::String,type::String,iz::Int64,particles::Vector{String},Ein::Vector{Float64},vec_Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,ΔEf::Float64,Ef⁻::Float64,Ef⁺::Float64,δi::Int64,Ui::Float64,Zi::Real,Ti::Float64,ri::Float64,subshells::Vector{String},Ec::Float64,gi::Int64,incoming_particle::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return dcs(interaction,L,Ei,Ef,type,gi,vec_Z,ωz,ρ,iz,Ein)
    elseif itype == Bremsstrahlung
        return dcs(interaction,L,Ei,Ef,Z,particle,type,iz,Ein,vec_Z,ωz,ρ)
    elseif itype == Compton
        return dcs(interaction,L,Ei,Ef,type,Ef⁻,Ef⁺,Z,iz)
    elseif itype == Elastic_Leptons
        return dcs(interaction,L,Ei,Z,particle,Ein[end],iz)
    elseif itype == Inelastic_Leptons
        return dcs(interaction,L,Ei,Ef,type,incoming_particle,Ui,Zi,Ti)
    elseif itype == Pair_Production
        return dcs(interaction,L,Ei,Ef,Z,particle,type,iz,particles,Ein[end])
    elseif itype == Photoelectric
        return dcs(interaction,L,Ei,Z,iz,δi,Ef⁻,Ef⁺)
    elseif itype == Rayleigh
        return dcs(interaction,L,Ei,Z,iz)
    elseif itype == Fluorescence
        return dcs(interaction,L,Ei,Z,iz,δi,Ef⁻,Ef⁺,incoming_particle)
    elseif itype == Auger
        return dcs(interaction,L,Ei,Z,iz,δi,Ef⁻,Ef⁺,incoming_particle)
    else
        error("Unknown interaction.")
    end
end

function tcs_dispatch(interaction::Interaction,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::String,Ecutoff::Float64,Eout::Vector{Float64},vec_Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,type::String)
    itype = typeof(interaction)
    if itype == Annihilation
        return tcs(interaction,Ei,Z)
    elseif itype == Bremsstrahlung
        return tcs(interaction,Ei,Z,Ec,iz,particle,type,Eout)
    elseif itype == Compton
        return tcs(interaction,Ei,Z,Eout,iz)
    elseif itype == Elastic_Leptons
        return tcs(interaction,Ei,Z,particle,Ecutoff,iz)
    elseif itype == Inelastic_Leptons
        return tcs(interaction,Ei,Ec,particle,Z)
    elseif itype == Pair_Production
        return tcs(interaction,Ei,Z,iz,Eout,type)
    elseif itype == Photoelectric
        return tcs(interaction,Ei,Z,iz)
    elseif itype == Rayleigh
        return tcs(interaction,Ei,Z,iz)
    else
        error("Unknown interaction.")
    end
end

function sp_dispatch(interaction::Interaction,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::String,Ecutoff::Float64,Eout::Vector{Float64})
    itype = typeof(interaction)
    if itype == Bremsstrahlung
        return sp(interaction,Z,ωz,ρ,state_of_matter,Ei,Ec,Eout)
    elseif itype == Inelastic_Leptons
        return sp(interaction,Z,ωz,ρ,state_of_matter,Ei,Ec,particle)
    else
        error("Unknown interaction.")
    end
end

function mt_dispatch(interaction::Interaction,Ei::Float64,Ec::Float64)
    itype = typeof(interaction)
    if itype == Bremsstrahlung
        return mt(interaction)
    elseif itype == Inelastic_Leptons
        return mt(interaction)
    else
        error("Unknown interaction.")
    end
end

function preload_data_dispatch(interaction::Interaction,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64,ωz::Vector{Float64},ρ::Float64,Eout::Vector{Float64},particle::String,type::String,Ein::Vector{Float64})
    itype = typeof(interaction)
    if itype == Annihilation
        return  preload_data(interaction,Z,Emax,Emin,L,type,Eout,Ein)
    elseif itype == Bremsstrahlung
        return preload_data(interaction,Z,Emax,Emin,L)
    elseif itype == Compton
        return preload_data(interaction,L,Z)
    elseif itype == Elastic_Leptons
        return preload_data(interaction,Z,L,particle)
    elseif itype == Inelastic_Leptons
        return preload_data(interaction,Z,ωz,ρ,particle)
    elseif itype == Pair_Production
        return preload_data(interaction,Z,Emax,Emin,Eout,L)
    elseif itype == Photoelectric
        return preload_data(interaction,Z,ωz,ρ,L)
    elseif itype == Rayleigh
        return preload_data(interaction,Z)
    elseif itype == Fluorescence
        return preload_data(interaction,Z,ωz,ρ,L,particle,Eout[end])
    elseif itype == Auger
        return preload_data(interaction,Z,ωz,ρ,L,particle,Eout[end])
    else
        error("Unknown interaction.")
    end
end

