"""
    particle_charge(particle::String)

Charge of the particle.

# Input Argument(s)
- 'particle::String': particle name.

# Output Argument(s)
- 'charge::Float64': particle charge.

# Author(s)
Charles Bienvenue

# Reference(s)

"""
function particle_charge(particle::String)
    charges = Dict("photons" => 0.0, "electrons" => -1.0, "positrons" => 1.0)
    if ~haskey(charges,particle) error("Unknown particle.") end
    return charges[particle]
end