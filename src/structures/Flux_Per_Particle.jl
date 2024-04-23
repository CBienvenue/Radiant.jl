
mutable struct Flux_Per_Particle

    # Variable(s)
    particle                        ::String
    flux                            ::Vector{Array{Float64,6}}
    flux_cutoff                     ::Vector{Array{Float64,5}}

    # Constructor(s)
    function Flux_Per_Particle(particle::String)

        this = new()
        this.particle = particle
        this.flux = Vector{Array{Float64,6}}()
        this.flux_cutoff = Vector{Array{Float64,5}}()

        return this
    end
end

# Method(s)
Base.propertynames(::Flux_Per_Particle) = 
(
    fieldnames(Flux_Per_Particle)...,
    :add_flux,
    :add_flux_cutoff,
    :get_flux,
    :get_flux_cutoff
)

function add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})
    push!(this.flux,flux)
end

function add_flux(this::Flux_Per_Particle,flux_per_particle::Flux_Per_Particle)
    if this.particle != flux_per_particle.particle error("Flux particle don't fit.") end
    append!(this.flux,flux_per_particle.flux)
    append!(this.flux_cutoff,flux_per_particle.flux_cutoff)
end

function add_flux_cutoff(this::Flux_Per_Particle,flux_cutoff::Array{Float64,5})
    push!(this.flux_cutoff,flux_cutoff)
end

function get_flux(this::Flux_Per_Particle)
    return sum(this.flux)
end

function get_flux_cutoff(this::Flux_Per_Particle)
    return sum(this.flux_cutoff)
end