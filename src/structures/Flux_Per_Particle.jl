"""
    Flux_Per_Particle

Structure used to contain flux information per particle.

"""
mutable struct Flux_Per_Particle

    # Variable(s)
    particle                        ::Particle
    flux                            ::Vector{Array{Float64,6}}
    flux_cutoff                     ::Vector{Array{Float64,5}}
    spectral_radius                 ::Vector{Vector{Float64}}

    # Constructor(s)
    function Flux_Per_Particle(particle::Particle)

        this = new()
        this.particle = particle
        this.flux = Vector{Array{Float64,6}}()
        this.flux_cutoff = Vector{Array{Float64,5}}()
        this.spectral_radius = Vector{Vector{Float64}}()

        return this
    end
end

# Method(s)
"""
    add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})

Add flux solution to the list of flux solutions.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux::Array{Float64,6}` : flux solution.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})
    push!(this.flux,flux)
end

"""
    add_flux(this::Flux_Per_Particle,flux_per_particle::Flux_Per_Particle)

Add flux solutions from another Flux_Per_Particle structure.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux_per_particle::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux_Per_Particle,flux_per_particle::Flux_Per_Particle)
    if get_tag(this.particle) != get_tag(flux_per_particle.particle) error("Flux particle don't fit.") end
    append!(this.flux,flux_per_particle.flux)
    append!(this.flux_cutoff,flux_per_particle.flux_cutoff)
    append!(this.spectral_radius,flux_per_particle.spectral_radius)
end

"""
    add_flux_cutoff(this::Flux_Per_Particle,flux_cutoff::Array{Float64,5})

Add flux at cutoff solutions to the list of flux at cutoff solutions.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `flux_cutoff::Array{Float64,5}` : flux at cutoff solution.

# Output Argument(s)
N/A

"""
function add_flux_cutoff(this::Flux_Per_Particle,flux_cutoff::Array{Float64,5})
    push!(this.flux_cutoff,flux_cutoff)
end

"""
    add_spectral_radius(this::Flux_Per_Particle,spectral_radius::Vector{Float64})

Add the per-energy-group in-group spectral-radius estimate of one generation to the list.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.
- `spectral_radius::Vector{Float64}` : estimated in-group spectral radius per energy group.

# Output Argument(s)
N/A

"""
function add_spectral_radius(this::Flux_Per_Particle,spectral_radius::Vector{Float64})
    push!(this.spectral_radius,spectral_radius)
end

"""
    get_spectral_radius(this::Flux_Per_Particle)

Get the per-energy-group in-group spectral-radius estimate of the last generation solved for the
particle (see the convergence-acceleration section of the documentation for its definition).

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `spectral_radius::Vector{Float64}` : estimated in-group spectral radius per energy group.

"""
function get_spectral_radius(this::Flux_Per_Particle)
    if length(this.spectral_radius) == 0 error("No spectral-radius data available.") end
    return this.spectral_radius[end]
end

"""
    get_flux(this::Flux_Per_Particle)

Get the total flux solution for the particle.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `flux::Array{Float64}` : flux solution.

"""
function get_flux(this::Flux_Per_Particle)
    return sum(this.flux)
end

"""
    get_flux_cutoff(this::Flux_Per_Particle)

Get the total flux solution at cutoff for the particle.

# Input Argument(s)
- `this::Flux_Per_Particle` : structure to contain flux solutions.

# Output Argument(s)
- `flux_cutoff::Array{Float64}` : flux solution at cutoff.

"""
function get_flux_cutoff(this::Flux_Per_Particle)
    return sum(this.flux_cutoff)
end