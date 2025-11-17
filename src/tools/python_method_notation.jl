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
    Particle,
    Spherical_Harmonics
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