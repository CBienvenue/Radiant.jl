"""
Cross_Sections

Structure used to define the parameters to extract or build a multigroup cross-sections library.

# Mandatory field(s)
- `source::String` : source of the cross-sections.
- `materials::Vector{Material}` : material list.
- **if `source = "fmac-m"`**
    - `file::String` : file containing cross-sections data.
- **if `source = "physics-models"`**
    - `particles::Vector{String}` : particle list.
    - `energy::Float64` : midpoint energy of the highest energy group [in MeV].
    - `number_of_groups::Int64` : number of energy groups.
    - `group_structure::String="log"` : type of group discretization.
    - `legendre_order::Int64` : maximum order of the angular Legendre moments of the differential cross-sections.
    - `interactions::Vector{Interaction}` : list of interaction.
- ** if `source = "custom"`**
    - `particles::Vector{String}` : particle list.


# Optional field(s) - with default values
- **if `source = "physics-models"`**
    - `cutoff::Float64 = 0.001` : lower energy bound of the lowest energy group (cutoff energy) [in MeV].
- ** if `source = "custom"`**
    - `custom_absorption::Vector{Real}` : absorption cross-sections per material.
    - `custom_scattering::Vector{Real}` : scattering (isotropic) cross-sections per material.

"""
mutable struct Cross_Sections

    # Variable(s)
    source                    ::Union{Missing,String}
    file                      ::Union{Nothing,String}
    number_of_materials       ::Int64
    materials                 ::Vector{Material}
    number_of_particles       ::Int64
    particles                 ::Vector{Particle}
    energy                    ::Union{Missing,Float64}
    cutoff                    ::Union{Missing,Float64}
    number_of_groups          ::Union{Missing,Vector{Int64}}
    group_structure           ::Union{Missing,Vector{String},String}
    interactions              ::Union{Missing,Vector{Interaction}}
    legendre_order            ::Union{Missing,Int64}
    energy_boundaries         ::Union{Missing,Vector{Vector{Float64}}}
    multigroup_cross_sections ::Union{Missing,Array{Multigroup_Cross_Sections}}
    is_build                  ::Bool
    custom_absorption         ::Vector{Real}
    custom_scattering         ::Vector{Real}

    # Constructor(s)
    function Cross_Sections()

        this = new()

        this.source = missing
        this.file = nothing
        this.number_of_materials = 0
        this.number_of_particles = 0
        this.materials = Vector{Material}()
        this.particles = Vector{String}()
        this.energy = missing
        this.cutoff = 0.001
        this.number_of_groups = missing
        this.group_structure = missing
        this.interactions = missing
        this.legendre_order = missing
        this.is_build = false
        this.multigroup_cross_sections = missing
        this.energy_boundaries = missing

        return this
    end
end

# Method(s)
"""
    is_ready_to_build(this::Cross_Sections)

Verify if all the required information to build the cross-section library is present.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
N/A

"""
function is_ready_to_build(this::Cross_Sections)
    if ismissing(this.source) error("Cannot build multigroup cross-sections data. The source of cross-sections data is not specified.") end
    if lowercase(this.source) == "fmac-m"
        if ismissing(this.file) error("Cannot build multigroup cross-sections data. The FMAC-M file name and its directory are not specified.") end
    elseif lowercase(this.source) == "physics-models"
        if ismissing(this.energy) error("Cannot build multigroup cross-sections data. The midpoint energy of the highest energy group is not specified.") end
        if ismissing(this.cutoff) error("Cannot build multigroup cross-sections data. The cutoff energy (lower bound of the lowest energy group) is not specified.") end
        if ismissing(this.number_of_groups) error("Cannot build multigroup cross-sections data. The number of energy group is not specified.") end
        if ismissing(this.group_structure) error("Cannot build multigroup cross-sections data. The multigroup structure is not specified.") end
        if ismissing(this.interactions) error("Cannot build multigroup cross-sections data. The type of interaction(s) between particle(s) and the material are not specified.") end
        if ismissing(this.legendre_order) error("Cannot build multigroup cross-sections data. The order for the Legendre expansion of differential scattering cross-sections is not specified.") end
        if length(this.number_of_groups) == 1
            this.number_of_groups = fill(this.number_of_groups[1],this.number_of_particles)
        elseif length(this.number_of_groups) != this.number_of_particles
            error("The number of groups information do not fit the number of particles.")
        end
        if length(this.group_structure) == 1
            this.group_structure = fill(this.group_structure[1],this.number_of_particles)
        elseif length(this.group_structure) != this.number_of_particles
            error("The number of groups structure information do not fit the number of particles.")
        end
    elseif lowercase(this.source) == "custom"
        ###
    else
        error("Unkown source of cross-sections data.")
    end
    if this.number_of_materials == 0 error("Cannot build multigroup cross-sections data. The material names are not specified.") end
    if this.number_of_particles == 0 error("Cannot build multigroup cross-sections data. The particle names are not specified.") end
end

"""
    build(this::Cross_Sections)

To build the cross-section library.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> ... # Defining the cross-sections library properties
julia> cs.build()
```
"""
function build(this::Cross_Sections)

    # Verification step
    is_ready_to_build(this)

    # Extract or produce formatted cross-sections data for transport calculations
    if lowercase(this.source) == "fmac-m"
        read_fmac_m(this)
    elseif lowercase(this.source) == "physics-models"
        try
            generate_cross_sections(this)
        catch
            rethrow()
        finally
            empty!(cache_radiant[]) # Clean cache before throwing an error.
        end
        empty!(cache_radiant[]) # Clean cache
    elseif lowercase(this.source) == "custom"
        custom_cross_sections(this)
    else
        error("Unkown source of cross-sections data.")
    end
    this.is_build = true

end

"""
    write(this::Cross_Sections,file::String)

To write a FMAC-M formatted file containing the cross-sections library.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `file::String` : file name and directory.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> ... # Defining the cross-sections library properties
julia> cs.write("fmac_m.txt")
```
"""
function write(this::Cross_Sections,file::String)
    if ~this.is_build this.build() end
    write_fmac_m(this,file)
end

"""
    set_source(this::Cross_Sections,source::String)

To define the source of the cross-sections library.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `source::String` : source of the cross-sections library which is either:
    - `source = "physics-models"` : multigroup cross-sections are produced by Radiant
    - `source = "FMAC-M"` : multigroup cross-sections are extracted from FMAC-M file.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_source("FMAC-M")
```
"""
function set_source(this::Cross_Sections,source::String)
    if lowercase(source) ∉ ["fmac-m","physics-models","custom"] error("Unkown source of cross-sections data.") end
    this.source = source
end

"""
    set_file(this::Cross_Sections,file::String)

To read a FMAC-M formatted file containing the cross-sections library.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `file::String` : file name and directory.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_file("fmac_m.txt")
```
"""
function set_file(this::Cross_Sections,file::String)
    this.file = file
end

"""
    set_materials(this::Cross_Sections,materials::Vector{Material})

To set the list of material, either contained in FMAC-M file in order, or to produce in Radiant.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `materials::Vector{Material}` : material list.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat1 = Material(); mat2 = Material()
julia> ... # Defining the material properties
julia> cs = Cross_Sections()
julia> cs.set_materials([mat1,mat2])
```
"""
function set_materials(this::Cross_Sections,materials::Union{Vector{Material},Material})
    if typeof(materials) == Material materials = [materials] end
    if length(materials) == 0 error("At least one material should be provided.") end
    this.materials = materials
    this.number_of_materials += length(materials)
end

"""
    set_particles(this::Cross_Sections,particles::Union{Vector{Particle},Particle})

To set the list of particles for which to produce coupled library of cross-sections.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particles::Vector{Particle}` : particles list.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles([electron,photon,positron])
```
"""
function set_particles(this::Cross_Sections,particles::Union{Vector{Particle},Particle})
    if typeof(particles) == Particle particles = [particles] end 
    if length(particles) == 0 error("At least one particle should be provided.") end
    this.particles = particles
    this.number_of_particles += length(particles)
end

"""
    set_energy(this::Cross_Sections,energy::Real)

To set the midpoint energy of the highest energy group.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `energy::Real` : midpoint energy of the highest energy group.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_energy(3.0)
```
"""
function set_energy(this::Cross_Sections,energy::Real)
    if energy < 0 error("The midpoint of the highest energy group has to be positive.") end
    this.energy = energy
end

"""
    set_cutoff(this::Cross_Sections,cutoff::Real)

To set the cutoff energy (lower bound of the lowest energy group).

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `cutoff::Real` : cutoff energy (lower bound of the lowest energy group)

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_cutoff(0.05)
```
"""
function set_cutoff(this::Cross_Sections,cutoff::Real)
    if cutoff < 0 error("The cutoff energy (lower bound of the lowest energy group) has to be positive.") end
    this.cutoff = cutoff
end

"""
    set_number_of_groups(this::Cross_Sections,number_of_groups::Vector{Int64})

To set the number of energy groups per particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `number_of_groups::Vector{Int64}` : number of energy groups per particle in order with the particle list.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles([electron,photon,positron])
julia> cs.set_number_of_groups([80,20,80]) # 80 groups with electrons and positrons, 20 with photons
```
"""
function set_number_of_groups(this::Cross_Sections,number_of_groups::Union{Vector{Int64},Int64})
    if (typeof(number_of_groups) == Int64) number_of_groups = [number_of_groups] end
    for g in number_of_groups if g ≤ 0 error("The number of energy groups should be at least one.") end end
    this.number_of_groups =  number_of_groups
end

"""
    set_group_structure(this::Cross_Sections,group_structure::Union{Vector{String},String})

To set the type of energy discretization structure per particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `group_structure::Vector{String}` : type of energy discretization structure per particle, where value per particle can take the following value:
    - `group_structure[i] = "linear"` : linearly spaced discretization.
    - `group_structure[i] = "log"` : logarithmically spaced discretization.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles([electron,photon,positron])
julia> cs.set_group_structure(["log","linear","log"]) # 80 groups with electrons and positrons, 20 with photons
```
"""
function set_group_structure(this::Cross_Sections,group_structure::Union{Vector{String},String})
    if typeof(group_structure) == String group_structure = [group_structure] end 
    for g in group_structure if g ∉ ["linear","log"] error("The group structure are either linearly or logarithmically spaced.") end end
    this.group_structure = group_structure
end

"""
    set_energy_boundaries(this::Cross_Sections,energy_boundaries::Vector{Vector{Float64}})

To set the energy boundaries per particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `energy_boundaries::Vector{Vector{Float64}}` : energy boundaries per particle.

# Output Argument(s)
N/A

"""
function set_energy_boundaries(this::Cross_Sections,energy_boundaries::Vector{Vector{Float64}})
    this.energy_boundaries = energy_boundaries
end

"""
    set_interactions(this::Cross_Sections,interactions::Vector{Interaction})

To set the interaction to take into account in the library of cross-sections.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `interactions::Vector{Interaction}` : list of interactions to use in the production of the cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles([electron])
julia> cs.set_interactions([Elastic_Collision(),Inelastic_Collision(),Bremsstrahlung(), Auger()])
```
"""
function set_interactions(this::Cross_Sections,interactions::Union{Vector{<:Interaction},<:Interaction})
    if (typeof(interactions) <: Interaction) interactions = [interactions] end
    if length(interactions) == 0 error("At least one interaction should be provided.") end
    this.interactions = interactions
end

"""
    set_legendre_order(this::Cross_Sections,legendre_order::Int64)

To set the maximum order of the Legendre expansion of the differential cross-sections.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `legendre_order::Int64` : maximum order of the Legendre expansion of the differential cross-sections.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_legendre_order(7)
```
"""
function set_legendre_order(this::Cross_Sections,legendre_order::Int64)
    if legendre_order < 0 error("Polynomial Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
end

"""
    set_multigroup_cross_sections(this::Cross_Sections,multigroup_cross_sections::Array{Multigroup_Cross_Sections})

To set the multigroup cross-sections libraries.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `multigroup_cross_sections::Array{Multigroup_Cross_Sections}` : multigroup cross-sections libraries.

# Output Argument(s)
N/A

"""
function set_multigroup_cross_sections(this::Cross_Sections,multigroup_cross_sections::Array{Multigroup_Cross_Sections})
    this.multigroup_cross_sections = multigroup_cross_sections
end

"""
    set_absorption(this::Cross_Sections,Σa::Vector{Float64})

To set absorption cross-sections.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `Σa::Vector{Float64}` : absorption cross-sections.

# Output Argument(s)
N/A

"""
function set_absorption(this::Cross_Sections,Σa::Vector{Float64})
    this.custom_absorption = Σa
end

"""
    set_scattering(this::Cross_Sections,Σs::Vector{Float64})

To set isotropic scattering cross-sections.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `Σs::Vector{Float64}` : scattering cross-sections.

# Output Argument(s)
N/A

"""
function set_scattering(this::Cross_Sections,Σs::Vector{Float64})
    this.custom_scattering = Σs
end

"""
    get_file(this::Cross_Sections)

Get the file directory and name.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `file::String` : file directory and name.

"""
function get_file(this::Cross_Sections)
    return this.file
end

"""
    get_number_of_materials(this::Cross_Sections)

Get the number of material.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `number_of_materials::Int64` : number of material.

"""
function get_number_of_materials(this::Cross_Sections)
    return this.number_of_materials
end

"""
    get_materials(this::Cross_Sections)

Get the material list.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `materials::Vector{Material}` : material list.

"""
function get_materials(this::Cross_Sections)
    return this.materials
end

"""
    get_number_of_particles(this::Cross_Sections)

Get the number of particles.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `number_of_particles::Int64` : number of particles.

"""
function get_number_of_particles(this::Cross_Sections)
    return this.number_of_particles
end

"""
    get_particles(this::Cross_Sections)

Get the particle list.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `particles::Vector{Particle}` : particle list.

"""
function get_particles(this::Cross_Sections)
    return this.particles
end

"""
    get_number_of_groups(this::Cross_Sections)

Get the number of groups for each particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `number_of_groups::Vector{Int64}` : number of groups for each particle.

"""
function get_number_of_groups(this::Cross_Sections)
    return this.number_of_groups
end

"""
    get_number_of_groups(this::Cross_Sections,particle::Particle)

Get the number of groups for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `number_of_groups::Int64` : number of groups.

"""
function get_number_of_groups(this::Cross_Sections,particle::Particle)
    index = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index) error("Cross-sections don't contain data for the given particle.") end
    return this.number_of_groups[index]
end

"""
    get_energy_boundaries(this::Cross_Sections,particle::Particle)

Get the energy boundaries for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `energy_boundaries::Vector{Float64}` : energy boundaries.

"""
function get_energy_boundaries(this::Cross_Sections,particle::Particle)
    index = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index) error("Cross-sections don't contain data for the given particle.") end
    return this.energy_boundaries[index]
end

"""
    get_energies(this::Cross_Sections,particle::Particle)

Get the group midpoint energies for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `energies::Vector{Float64}` : group midpoint energies.

"""
function get_energies(this::Cross_Sections,particle::Particle)
    Eb = this.get_energy_boundaries(particle)
    return (Eb[1:end-1] + Eb[2:end])/2
end

"""
    get_energy_width(this::Cross_Sections,particle::Particle)

Get the energy group width for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `width::Vector{Float64}` : energy group width.

"""
function get_energy_width(this::Cross_Sections,particle::Particle)
    Eb = this.get_energy_boundaries(particle)
    return (Eb[1:end-1] - Eb[2:end])
end

"""
    get_absorption(this::Cross_Sections,particle::Particle)

Get the absorption cross-sections for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `Σa::Array{Float64}` : absorption cross-sections.

"""
function get_absorption(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup cross-sections. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σa = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σa[:,n] = this.multigroup_cross_sections[index_particle,n].get_absorption()
    end
    return Σa
end

"""
    get_total(this::Cross_Sections,particle::Particle)

Get the total cross-sections for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `Σt::Array{Float64}` : total cross-sections.

"""
function get_total(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup cross-sections. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σt = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σt[:,n] = this.multigroup_cross_sections[index_particle,n].get_total()
    end
    return Σt
end

"""
    get_scattering(this::Cross_Sections,particle_in::Particle,particle_out::Particle,
    legendre_order::Int64)

Get the Legendre moments of the scattering cross-sections for a specified incoming and
outgoing particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle_in::Particle` : incoming particle.
- `particle_out::Particle` : outgoing particle.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Σs::Array{Float64}` : Legendre moments of the scattering cross-sections.

"""
function get_scattering(this::Cross_Sections,particle_in::Particle,particle_out::Particle,legendre_order::Int64)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup scattering cross-sections. Missing data.") end
    index_particle_in = findfirst(x -> x == particle_in,this.get_particles())
    index_particle_out = findfirst(x -> x == particle_out,this.get_particles())
    if isnothing(index_particle_in) || isnothing(index_particle_out) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ngi = this.get_number_of_groups(particle_in)
    Ngf = this.get_number_of_groups(particle_out)
    Σs = zeros(Nmat,Ngi,Ngf,legendre_order+1)
    for n in range(1,Nmat)
        Σs_i = this.multigroup_cross_sections[index_particle_in,n].get_scattering(index_particle_out)
        Σs[n,:,:,1:min(legendre_order+1,length(Σs_i[1,1,:]))] = Σs_i[:,:,1:min(legendre_order+1,length(Σs_i[1,1,:]))]
    end
    return Σs
end

"""
    get_boundary_stopping_powers(this::Cross_Sections,particle::Particle)

Get the stopping powers at group boundaries for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `S::Array{Float64}` : stopping powers at group boundaries.

"""
function get_boundary_stopping_powers(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup stopping powers. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Sb = zeros(Ng+1,Nmat)
    for n in range(1,Nmat)
        Sb[:,n] = this.multigroup_cross_sections[index_particle,n].get_boundary_stopping_powers()
    end
    return Sb
end

"""
    get_stopping_powers(this::Cross_Sections,particle::Particle)

Get the stopping powers for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `S::Array{Float64}` : stopping powers.

"""
function get_stopping_powers(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup stopping powers. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    S = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        S[:,n] = this.multigroup_cross_sections[index_particle,n].get_stopping_powers()
    end
    return S
end

"""
    get_momentum_transfer(this::Cross_Sections,particle::Particle)

Get the momentum transfers for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `T::Array{Float64}` : momentum transfers.

"""
function get_momentum_transfer(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    T = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        T[:,n] = this.multigroup_cross_sections[index_particle,n].get_momentum_transfer()
    end
    return T
end

"""
    get_energy_deposition(this::Cross_Sections,particle::Particle)

Get the energy deposition cross-sections for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `Σe::Array{Float64}` : energy deposition cross-sections.

"""
function get_energy_deposition(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σe = zeros(Ng+1,Nmat)
    for n in range(1,Nmat)
        Σe[:,n] = this.multigroup_cross_sections[index_particle,n].get_energy_deposition()
    end
    return Σe
end

"""
    get_charge_deposition(this::Cross_Sections,particle::Particle)

Get the charge deposition cross-sections for a specified particle.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.
- `particle::Particle` : particle.

# Output Argument(s)
- `Σe::Array{Float64}` : charge deposition cross-sections.

"""
function get_charge_deposition(this::Cross_Sections,particle::Particle)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    if isnothing(index_particle) error("Cross-sections don't contain data for the given particle.") end
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σc = zeros(Ng+1,Nmat)
    for n in range(1,Nmat)
        Σc[:,n] = this.multigroup_cross_sections[index_particle,n].get_charge_deposition()
    end
    return Σc
end

"""
    get_densities(this::Cross_Sections)

Get the material densities.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `ρ::Vector{Float64}` : densities.

"""
function get_densities(this::Cross_Sections)
    Nmat = this.get_number_of_materials()
    ρ = zeros(Nmat)
    for n in range(1,Nmat)
        ρ[n] = this.materials[n].get_density()
    end
    return ρ
end

"""
    get_energy(this::Cross_Sections)

Get the midpoint energy of the highest energy group.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `energy::Float64` : midpoint energy of the highest energy group.

"""
function get_energy(this::Cross_Sections)
    return this.energy
end

"""
    get_cutoff(this::Cross_Sections)

Get the cutoff energy.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `cutoff::Float64` : cutoff energy.

"""
function get_cutoff(this::Cross_Sections)
    return this.cutoff
end

"""
    get_legendre_order(this::Cross_Sections)

Get the Legendre truncation order.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `legendre_order::Int64` : Legendre truncation order.

"""
function get_legendre_order(this::Cross_Sections)
    return this.legendre_order
end

"""
    get_group_structure(this::Cross_Sections)

Get the type of group structure.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `group_structure::String` : type of group structure.

"""
function get_group_structure(this::Cross_Sections)
    return this.group_structure
end

"""
    get_interactions(this::Cross_Sections)

Get the interaction list.

# Input Argument(s)
- `this::Cross_Sections` : cross-sections library.

# Output Argument(s)
- `interactions::Vector{Interation}` : interaction list.

"""
function get_interactions(this::Cross_Sections)
    return this.interactions
end