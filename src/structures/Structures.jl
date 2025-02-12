# Data structures

include("Particle.jl")
export Particle, Photon, Electron, Positron

# To stock information about cross section generation
include("Interaction.jl")
include("Elastic_Leptons.jl")
include("Inelastic_Leptons.jl")
include("Bremsstrahlung.jl")
include("Compton.jl")
include("Pair_Production.jl")
include("Photoelectric.jl")
include("Annihilation.jl")
include("Rayleigh.jl")
include("Fluorescence.jl")
include("Auger.jl")
export Elastic_Leptons
export Inelastic_Leptons
export Bremsstrahlung
export Compton
export Pair_Production
export Photoelectric
export Annihilation
export Rayleigh
export Fluorescence
export Auger

# To organize information for transport calculations
include("Multigroup_Cross_Sections.jl")
include("Material.jl")
include("Cross_Sections.jl")
include("Geometry.jl")
include("Discrete_Ordinates.jl")
include("Solvers.jl")
include("Surface_Source.jl")
include("Volume_Source.jl")
include("Source.jl")
include("Sources.jl")
include("Fixed_Sources.jl")
include("Flux_Per_Particle.jl")
include("Flux.jl")
include("Computation_Unit.jl")
export Material
export Cross_Sections
export Geometry
export Discrete_Ordinates
export Solvers
export Surface_Source
export Volume_Source
export Fixed_Sources
export Computation_Unit

# To use Python-like notation for methods
RadiantObject = Union{
    Material,
    Cross_Sections,
    Multigroup_Cross_Sections,
    Geometry,
    Discrete_Ordinates,
    Solvers,
    Surface_Source,
    Volume_Source,
    Fixed_Sources,
    Source,
    Sources,
    Computation_Unit,
    Flux_Per_Particle,
    Flux,
    Interaction,
    Particle
}

"""
    Base.getproperty(object::RadiantObject, s::Symbol)

Special function to enable calling a method "m(this::T,...)" on a struct "S" of type "T"
using the notation "S.m(...)" for a subset of struct.

"""
function Base.getproperty(object::RadiantObject, s::Symbol)
    if hasfield(typeof(object), s)
        return getfield(object, s)
    else
        # Check if the function is defined in the Radiant module
        if isdefined(Radiant, s) && getproperty(Radiant, s) isa Function
            try
                return function(x...)
                    func = getproperty(Radiant, s)
                    return func(object, x...)
                end
            catch e
                # Log the original error and rethrow it
                type = typeof(object)
                @error "Failed to invoke function '$s' on object of type '$type'. Error: $e"
                rethrow()
            end
        else
            type = typeof(object)
            throw(error("The object of type '$type' has no function '$s'."))
        end
    end
end