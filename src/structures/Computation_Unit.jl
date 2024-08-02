"""
    Computation_Unit

Structure used to consolidate the cross-sections, geometry, solvers and sources, execute transport calculations and extract its results.

# Mandatory field(s)
- `cross_sections::Cross_Sections`: cross-section library.
- `geometry::Geometry`: geometry.
- `solvers::Solvers`: solvers.
- `sources::Sources`: fixed sources.

# Optional field(s) - with default values
- N/A

"""
mutable struct Computation_Unit

    # Variable(s)
    cross_sections            ::Union{Missing,Cross_Sections}
    geometry                  ::Union{Missing,Geometry}
    solvers                   ::Union{Missing,Solvers}
    sources                   ::Union{Missing,Fixed_Sources}
    flux                      ::Union{Missing,Flux}

    # Constructor(s)
    function Computation_Unit()

        this = new()

        this.cross_sections = missing
        this.geometry = missing
        this.solvers = missing
        this.sources = missing
        this.flux = missing

        return this
    end
end

# Method(s)
Base.propertynames(::Computation_Unit) = 
(
    fieldnames(Computation_Unit)...,
    :set_cross_sections,
    :set_geometry,
    :set_methods,
    :set_sources,
    :run,
    :get_flux,
    :get_energy_deposition,
    :get_charge_deposition,
    :get_voxels_position,
    :get_energies
)

"""
    set_cross_sections(this::Computation_Unit,cross_sections::Cross_Sections)

To set the cross-sections library for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `cross_sections::Cross_Sections`: cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> ... # Define cross-sections properties and generate multigroup cross-sections.
julia> cu = Computation_Unit()
julia> cu.set_cross_sections(cs)
```
"""
function set_cross_sections(this::Computation_Unit,cross_sections::Cross_Sections)
    this.cross_sections = cross_sections
end

"""
    set_geometry(this::Computation_Unit,geometry::Geometry)

To set the geometry for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `geometry::Geometry`: geometry.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> ... # Define geometry and its properties
julia> cu = Computation_Unit()
julia> cu.set_geometry(geo)
```
"""
function set_geometry(this::Computation_Unit,geometry::Geometry)
    this.geometry = geometry
end

"""
set_solvers(this::Computation_Unit,solvers::Solvers)

To set the solvers for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `solvers::Solvers`: collection of solvers per particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Solvers()
julia> ... # Define all the discretization solvers and their properties
julia> cu = Computation_Unit()
julia> cu.set_solvers(ms)
```
"""
function set_solvers(this::Computation_Unit,solvers::Solvers)
    this.solvers = solvers
end

"""
    set_sources(this::Computation_Unit,sources::Fixed_Sources)

To set the fixed sources for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `sources::Fixed_Sources`: collection of fixed sources.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> fs = Fixed_Sources()
julia> ... # Define all the fixed sources and their properties
julia> cu = Computation_Unit()
julia> cu.set_sources(fs)
```
"""
function set_sources(this::Computation_Unit,sources::Fixed_Sources)
    this.sources = sources
end

"""
    run(this::Computation_Unit)

To lauch the transport calculations to solve the transport equation and obtain the flux solution.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define the cross-sections, geometry, fixed sources and discretization solvers
julia> cu.run()
```
"""
function run(this::Computation_Unit)
    is_CUDA = false
    #reset_timer!()
    this.flux = transport(this.cross_sections,this.geometry,this.solvers,this.sources,is_CUDA)
    #print_timer()
end

"""
    get_energy_deposition(this::Computation_Unit,type::String)

To get the array containing the energy deposition in each voxels.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `type::String`: type of energy deposition, which can takes the following values:
    - `type = "total"`: total energy deposition.
    - `type = "photons"`: photons energy deposition.
    - `type = "electrons"`: electrons energy deposition.
    - `type = "positrons"`: positrons energy deposition.

# Output Argument(s)
- `energy_deposition::Array{Float64}`: energy deposition array.

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define computation unit and run it.
julia> energy_deposition = cu.get_energy_deposition("total")
```
"""
function get_energy_deposition(this::Computation_Unit,type::String)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract energy deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of energy deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Energy deposition for the specified particle is not available.") end
    return energy_deposition(this.cross_sections,this.geometry,this.solvers,this.sources,this.flux,type)
end

"""
    get_charge_deposition(this::Computation_Unit,type::String)

To get the array containing the charge deposition in each voxels.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `type::String`: type of charge deposition, which can takes the following values:
    - `type = "total"`: total charge deposition.
    - `type = "photons"`: photons charge deposition.
    - `type = "electrons"`: electrons charge deposition.
    - `type = "positrons"`: positrons charge deposition.

# Output Argument(s)
- `charge_deposition::Array{Float64}`: charge deposition array.

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define computation unit and run it.
julia> charge_deposition = cu.get_charge_deposition("total")
```
"""
function get_charge_deposition(this::Computation_Unit,type::String)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract charge deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of charge deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Charge deposition for the specified particle is not available.") end
    return charge_deposition(this.cross_sections,this.geometry,this.solvers,this.sources,this.flux,type)
end

"""
    get_flux(this::Computation_Unit,particle::String)

To get the array containing the flux in each voxels and in each energy group for the specified particle.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `particle::String`: flux particle, which can takes the following values:
    - `particle = "photons"`: photons charge deposition.
    - `particle = "electrons"`: electrons charge deposition.
    - `particle = "positrons"`: positrons charge deposition.

# Output Argument(s)
- `flux::Array{Float64}`: flux array.

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define computation unit and run it.
julia> flux = cu.get_flux("electrons")
```
"""
function get_flux(this::Computation_Unit,particle::String)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract flux, please use .run() method before.") end
    if particle ∉ ["electrons","photons","positrons"] error("Unknown particle for flux.") end
    if particle ∈ ["electrons","photons","positrons"] && particle ∉ this.flux.get_particles() error("Flux for the specified particle is not available.") end
    return flux(this.cross_sections,this.geometry,this.flux,particle)
end

"""
    get_voxels_position(this::Computation_Unit,axis::String)

To set the mid-point voxels position along the specified axis.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `axis::String`: axis, which can takes the following values:
    - `boundary = "x"`: along x-axis
    - `boundary = "y"`: along y-axis
    - `boundary = "z"`: along z-axis

# Output Argument(s)
- `x::Vector{Float64}`: mid-point voxels position along the specified axis.

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define computation unit and run it.
julia> x = cu.get_voxels_position("x")
```
"""
function get_voxels_position(this::Computation_Unit,axis::String)
    return this.geometry.get_voxels_position(axis)
end

"""
    get_energies(this::Computation_Unit,particle::String)

To set the mid-point energy in each group for the specified particle.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `particle::String`: particle identifier, where each particle is either:
    - `particle = "photons"`: photons.
    - `particle = "electrons"`: electrons.
    - `particle = "positrons"`: positrons.

# Output Argument(s)
- `E::Vector{Float64}`: mid-point energy in each group for the specified particle.

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define computation unit and run it.
julia> E = cu.get_voxels_position("electrons")
```
"""
function get_energies(this::Computation_Unit,particle::String)
    return this.cross_sections.get_energies(particle)
end