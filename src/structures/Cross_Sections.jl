"""
Cross_Sections

Structure used to define the parameters to extract or build a multigroup cross-sections library.

# Mandatory field(s)
- `name::String`: name (or identifier) of the Cross_Sections structure.
- `source::String`: source of the cross-sections.
- **if `source = "FMAC-M"`**
    - `file::String`: file containing cross-sections data.
    - `materials::Vector{Material}`: material list.
- **if `source = "RADIANT"`**
    - `materials::Vector{Material}`: material list.
    - `particles::Vector{String}`: particle list.
    - `energy::Float64`: midpoint energy of the highest energy group [in MeV].
    - `number_of_groups::Int64`: number of energy groups.
    - `group_structure::String="log"`: type of group discretization.
    - `legendre_order::Int64`: maximum order of the angular Legendre moments of the differential cross-sections.
    - `interactions::Vector{Interaction}`: list of interaction.
    - `solvers::Vector{String}`: list of solver per particle.

# Optional field(s) - with default values
- `cutoff::Float64=0.001`: lower energy bound of the lowest energy group (cutoff energy) [in MeV].

"""
mutable struct Cross_Sections

    # Variable(s)
    name                      ::Union{Missing,String}
    source                    ::Union{Missing,String}
    file                      ::Union{Nothing,String}
    number_of_materials       ::Int64
    materials                 ::Vector{Material}
    number_of_particles       ::Int64
    particles                 ::Vector{String}
    energy                    ::Union{Missing,Float64}
    cutoff                    ::Union{Missing,Float64}
    number_of_groups          ::Union{Missing,Vector{Int64}}
    group_structure           ::Union{Missing,Vector{String},String}
    interactions              ::Union{Missing,Vector{Interaction}}
    solvers                   ::Union{Missing,Vector{String},String}
    legendre_order            ::Union{Missing,Int64}
    energy_boundaries         ::Union{Missing,Vector{Vector{Float64}}}
    multigroup_cross_sections ::Union{Missing,Array{Multigroup_Cross_Sections}}
    is_build                  ::Bool

    # Constructor(s)
    function Cross_Sections()

        this = new()

        this.name = missing
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
        this.solvers = missing
        this.legendre_order = missing
        this.is_build = false
        this.multigroup_cross_sections = missing
        this.energy_boundaries = missing

        return this
    end
end

# Method(s)
Base.propertynames(::Cross_Sections) = 
(
    fieldnames(Cross_Sections)...,
    :set_source,
    :set_file,
    :set_materials,
    :set_particles,
    :set_energy,
    :set_cutoff,
    :set_number_of_groups,
    :set_group_structure,
    :set_energy_boundaries,
    :set_interactions,
    :set_solvers,
    :set_legendre_order,
    :set_multigroup_cross_sections,
    :get_file,
    :get_materials,
    :get_number_of_materials,
    :get_particles,
    :get_number_of_groups,
    :get_energy_boundaries,
    :get_energies,
    :get_energy_width,
    :get_total,
    :get_absorption,
    :get_scattering,
    :get_stopping_powers,
    :get_momentum_transfer,
    :get_energy_deposition,
    :get_charge_deposition,
    :get_densities,
    :get_energy,
    :get_cutoff,
    :get_legendre_order,
    :get_solvers,
    :get_group_structure,
    :get_interactions,
    :build
)

function is_ready_to_build(this::Cross_Sections)
    if ismissing(this.source) error("Cannot build multigroup cross-sections data. The source of cross-sections data is not specified.") end
    if uppercase(this.source) == "FMAC-M"
        if ismissing(this.file) error("Cannot build multigroup cross-sections data. The FMAC-M file name and its directory are not specified.") end
    elseif uppercase(this.source) == "RADIANT"
        if ismissing(this.energy) error("Cannot build multigroup cross-sections data. The midpoint energy of the highest energy group is not specified.") end
        if ismissing(this.cutoff) error("Cannot build multigroup cross-sections data. The cutoff energy (lower bound of the lowest energy group) is not specified.") end
        if this.number_of_particles == 0 error("Cannot build multigroup cross-sections data. The particle names are not specified.") end
        if ismissing(this.number_of_groups) error("Cannot build multigroup cross-sections data. The number of energy group is not specified.") end
        if ismissing(this.group_structure) error("Cannot build multigroup cross-sections data. The multigroup structure is not specified.") end
        if ismissing(this.interactions) error("Cannot build multigroup cross-sections data. The type of interaction(s) between particle(s) and the material are not specified.") end
        if ismissing(this.solvers) error("Cannot build multigroup cross-sections data. The type of solver(s) for each particle transport are not specified.") end
        if ismissing(this.legendre_order) error("Cannot build multigroup cross-sections data. The order for the Legendre expansion of differential scattering cross-sections is not specified.") end
    else
        error("Unkown source of cross-sections data.")
    end
    if this.number_of_materials == 0 error("Cannot build multigroup cross-sections data. The material names are not specified.") end
end

"""
    build(this::Cross_Sections)

To build the cross-section library.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.

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
    if uppercase(this.source) == "FMAC-M"
        read_fmac_m(this)
    elseif uppercase(this.source) == "RADIANT"
        generate_cross_sections(this)
    else
        error("Unkown source of cross-sections data.")
    end
    this.is_build = true

end

"""
    write(this::Cross_Sections,file::String)

To write a FMAC-M formatted file containing the cross-sections library.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `file::String`: file name and directory.

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
- `this::Cross_Sections`: cross-sections library.
- `source::String`: source of the cross-sections library which is either:
    - `source = "RADIANT"`: multigroup cross-sections are produced by Radiant
    - `source = "FMAC-M"`: multigroup cross-sections are extracted from FMAC-M file.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_source("FMAC-M")
```
"""
function set_source(this::Cross_Sections,source::String)
    if uppercase(source) ∉ ["FMAC-M","RADIANT"] error("Unkown source of cross-sections data.") end
    this.source = source
end

"""
    set_file(this::Cross_Sections,file::String)

To read a FMAC-M formatted file containing the cross-sections library.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `file::String`: file name and directory.

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
- `this::Cross_Sections`: cross-sections library.
- `materials::Vector{Material}`: material list.

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
function set_materials(this::Cross_Sections,materials::Vector{Material})
    if length(materials) == 0 error("At least one material should be provided.") end
    this.materials = materials
    this.number_of_materials += length(materials)
end

"""
    set_particles(this::Cross_Sections,particles::Vector{String})

To set the list of particles for which to produce coupled library of cross-sections.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `particles::Vector{String}`: particles list, where each particle is either:
    - `particles[i] = "photons"`: photons production and interaction are taken into account.
    - `particles[i] = "electrons"`: electrons production and interaction are taken into account.
    - `particles[i] = "positrons"`: positrons production and interaction are taken into account.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles(["electrons","photons","positrons"])
```
"""
function set_particles(this::Cross_Sections,particles::Vector{String})
    if length(particles) == 0 error("At least one particle should be provided.") end
    for p in particles if lowercase(p) ∉ ["photons","electrons","positrons"] error("Unknown particle type") end end
    this.particles = particles
    this.number_of_particles += length(particles)
end

"""
    set_energy(this::Cross_Sections,energy::Real)

To set the midpoint energy of the highest energy group.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `energy::Real`: midpoint energy of the highest energy group.

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
- `this::Cross_Sections`: cross-sections library.
- `cutoff::Real`: cutoff energy (lower bound of the lowest energy group)

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
- `this::Cross_Sections`: cross-sections library.
- `number_of_groups::Vector{Int64}`: number of energy groups per particle in order with the particle list.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles(["electrons","photons","positrons"])
julia> cs.set_number_of_groups([80,20,80]) # 80 groups with leptons, 20 with photons
```
"""
function set_number_of_groups(this::Cross_Sections,number_of_groups::Vector{Int64})
    for g in number_of_groups if g ≤ 0 error("The number of energy groups should be at least one.") end end
    this.number_of_groups = number_of_groups
end

"""
    set_group_structure(this::Cross_Sections,group_structure::Vector{String})

To set the type of energy discretization structure per particle.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `group_structure::Vector{String}`: type of energy discretization structure per particle, where value per particle can take the following value:
    - `group_structure[i] = "linear"`: linearly spaced discretization.
    - `group_structure[i] = "log"`: logarithmically spaced discretization.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles(["electrons","photons","positrons"])
julia> cs.set_group_structure(["log","linear","log"]) # 80 groups with leptons, 20 with photons
```
"""
function set_group_structure(this::Cross_Sections,group_structure::Vector{String})
    for g in group_structure if g ∉ ["linear","log"] error("The group structure are either linearly or logarithmically spaced.") end end
    this.group_structure = group_structure
end

function set_energy_boundaries(this::Cross_Sections,energy_boundaries::Vector{Vector{Float64}})
    this.energy_boundaries = energy_boundaries
end

"""
    set_interactions(this::Cross_Sections,interactions::Vector{Interaction})

To set the interaction to take into account in the library of cross-sections.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `interactions::Vector{Interaction}`: list of interactions to use in the production of the cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles(["electrons"])
julia> cs.set_interactions([Elastic_Leptons(),Inelastic_Leptons(),Bremsstrahlung(), Auger()])
```
"""
function set_interactions(this::Cross_Sections,interactions::Vector{Interaction})
    if length(interactions) == 0 error("At least one interaction should be provided.") end
    this.interactions = interactions
end

"""
    set_solvers(this::Cross_Sections,solvers::Vector{String})

To set the transport solver per particle.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `solvers::Vector{String}`: list of the transport to use per particle, which is defined by
    - `solvers[i] = "BTE"`: Boltzmann transport equation
    - `solvers[i] = "BFP"`: Boltzmann Fokker-Planck equation
    - `solvers[i] = "BCSD"`: Boltzmann-CSD equation
    - `solvers[i] = "FP"`: Fokker-Planck equation
    - `solvers[i] = "CSD"`: Continuous slowing-down only equation
    - `solvers[i] = "BFP-EF"`: Boltzmann Fokker-Planck without elastic

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> cs.set_particles(["electrons","photons"])
julia> cs.set_group_structure("BFP","BTE"]) # BFP for electrons, BTE for photons 
```
"""
function set_solvers(this::Cross_Sections,solvers::Vector{String})
    for solver in solvers if uppercase(solver) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("The solver type is either BTE, BFP, BCSD, FP, CSD or BFP-EF.") end end
    this.solvers = solvers
end

"""
    set_legendre_order(this::Cross_Sections,legendre_order::Int64)

To set the maximum order of the Legendre expansion of the differential cross-sections.

# Input Argument(s)
- `this::Cross_Sections`: cross-sections library.
- `legendre_order::Int64`: maximum order of the Legendre expansion of the differential cross-sections.

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

function set_multigroup_cross_sections(this::Cross_Sections,multigroup_cross_sections::Array{Multigroup_Cross_Sections})
    this.multigroup_cross_sections = multigroup_cross_sections
end

function get_file(this::Cross_Sections)
    return this.file
end

function get_number_of_materials(this::Cross_Sections)
    return this.number_of_materials
end

function get_materials(this::Cross_Sections)
    return this.materials
end

function get_number_of_particles(this::Cross_Sections)
    return this.number_of_particles
end

function get_particles(this::Cross_Sections)
    return this.particles
end

function get_number_of_groups(this::Cross_Sections)
    return this.number_of_groups
end

function get_number_of_groups(this::Cross_Sections,particle::String)
    index = findfirst(x -> x == particle,this.get_particles())
    return this.number_of_groups[index]
end

function get_energy_boundaries(this::Cross_Sections,particle::String)
    index = findfirst(x -> x == particle,this.get_particles())
    return this.energy_boundaries[index]
end

function get_energies(this::Cross_Sections,particle::String)
    Eb = this.get_energy_boundaries(particle)
    return (Eb[1:end-1] + Eb[2:end])/2
end

function get_energy_width(this::Cross_Sections,particle::String)
    Eb = this.get_energy_boundaries(particle)
    return (Eb[1:end-1] - Eb[2:end])
end

function get_absorption(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup cross-sections. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σa = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σa[:,n] = this.multigroup_cross_sections[index_particle,n].get_absorption()
    end
    return Σa
end

function get_total(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup cross-sections. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σt = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σt[:,n] = this.multigroup_cross_sections[index_particle,n].get_total()
    end
    return Σt
end

function get_scattering(this::Cross_Sections,particle_in::String,particle_out::String,legendre_order::Int64)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup scattering cross-sections. Missing data.") end
    index_particle_in = findfirst(x -> x == particle_in,this.get_particles())
    index_particle_out = findfirst(x -> x == particle_out,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle_in)
    Σs = zeros(Nmat,Ng,Ng,legendre_order+1)
    for n in range(1,Nmat)
        Σs_i = this.multigroup_cross_sections[index_particle_in,n].get_scattering(index_particle_out)
        if length(Σs_i[1,1,:]) < legendre_order error("Scattering cross-sections is not expanded up to ",legendre_order,"th-order Legendre moment.") end
        Σs[n,:,:,:] = Σs_i[:,:,1:legendre_order+1]
    end
    return Σs
end

function get_stopping_powers(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup stopping powers. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    β = zeros(Ng+1,Nmat)
    for n in range(1,Nmat)
        β[:,n] = this.multigroup_cross_sections[index_particle,n].get_stopping_powers()
    end
    return β
end

function get_momentum_transfer(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    α = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        α[:,n] = this.multigroup_cross_sections[index_particle,n].get_momentum_transfer()
    end
    return α
end

function get_energy_deposition(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σe = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σe[:,n] = this.multigroup_cross_sections[index_particle,n].get_energy_deposition()
    end
    return Σe
end

function get_charge_deposition(this::Cross_Sections,particle::String)
    if ismissing(this.multigroup_cross_sections) error("Unable to get multigroup momentum transfer. Missing data.") end
    index_particle = findfirst(x -> x == particle,this.get_particles())
    Nmat = this.get_number_of_materials()
    Ng = this.get_number_of_groups(particle)
    Σc = zeros(Ng,Nmat)
    for n in range(1,Nmat)
        Σc[:,n] = this.multigroup_cross_sections[index_particle,n].get_charge_deposition()
    end
    return Σc
end

function get_densities(this::Cross_Sections)
    Nmat = this.get_number_of_materials()
    ρ = zeros(Nmat)
    for n in range(1,Nmat)
        ρ[n] = this.materials[n].get_density()
    end
    return ρ
end

function get_energy(this::Cross_Sections)
    return this.energy
end

function get_cutoff(this::Cross_Sections)
    return this.cutoff
end

function get_legendre_order(this::Cross_Sections)
    return this.legendre_order
end

function get_solvers(this::Cross_Sections)
    return this.solvers
end

function get_group_structure(this::Cross_Sections)
    return this.group_structure
end

function get_interactions(this::Cross_Sections)
    return this.interactions
end