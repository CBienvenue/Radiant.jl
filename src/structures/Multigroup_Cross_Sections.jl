"""
    Multigroup_Cross_Sections

Structure used to multigroup cross-sections for a given particle.

"""
mutable struct Multigroup_Cross_Sections

    # Variable(s)
    number_of_groups               ::Int64
    total                          ::Union{Missing,Vector{Float64}}
    absorption                     ::Union{Missing,Vector{Float64}}
    boundary_stopping_powers       ::Union{Missing,Vector{Float64}}
    stopping_powers                ::Union{Missing,Vector{Float64}}
    momentum_transfer              ::Union{Missing,Vector{Float64}}
    energy_deposition              ::Union{Missing,Vector{Float64}}
    charge_deposition              ::Union{Missing,Vector{Float64}}
    scattering                     ::Vector{Array{Float64,3}}

    # Constructor(s)
    function Multigroup_Cross_Sections(number_of_groups::Int64)

        this = new()
        this.number_of_groups = number_of_groups
        this.total = missing
        this.absorption = missing
        this.boundary_stopping_powers = missing
        this.stopping_powers = missing
        this.momentum_transfer = missing
        this.energy_deposition = missing
        this.charge_deposition = missing
        this.scattering = Vector{Array{Float64,3}}()

        return this
    end
end

# Method(s)
"""
    set_total(this::Multigroup_Cross_Sections,total::Vector{Float64})

To set the total cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `total::Vector{Float64}` : total cross-sections.

# Output Argument(s)
N/A

"""
function set_total(this::Multigroup_Cross_Sections,total::Vector{Float64})
    if length(total) != this.number_of_groups error("The length of the total cross-sections don't fit the number of groups.") end
    this.total = total
end

"""
    set_absorption(this::Multigroup_Cross_Sections,absorption::Vector{Float64})

To set the absorption cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `absorption::Vector{Float64}` : absorption cross-sections.

# Output Argument(s)
N/A

"""
function set_absorption(this::Multigroup_Cross_Sections,absorption::Vector{Float64})
    if length(absorption) != this.number_of_groups error("The length of the absorption cross-sections don't fit the number of groups.") end
    this.absorption = absorption
end

"""
    set_boundary_stopping_powers(this::Multigroup_Cross_Sections,stopping_powers::Vector{Float64})

To set the stopping powers at group boundaries.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `boundary_stopping_powers::Vector{Float64}` : stopping powers at group boundaries.

# Output Argument(s)
N/A

"""
function set_boundary_stopping_powers(this::Multigroup_Cross_Sections,boundary_stopping_powers::Vector{Float64})
    if length(boundary_stopping_powers) != this.number_of_groups + 1 error("The length of the boundary_stopping_powers don't fit the number of groups.") end
    this.boundary_stopping_powers = boundary_stopping_powers
end

"""
    set_stopping_powers(this::Multigroup_Cross_Sections,stopping_powers::Vector{Float64})

To set the stopping powers.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `stopping_powers::Vector{Float64}` : stopping powers.

# Output Argument(s)
N/A

"""
function set_stopping_powers(this::Multigroup_Cross_Sections,stopping_powers::Vector{Float64})
    if length(stopping_powers) != this.number_of_groups error("The length of the stopping_powers don't fit the number of groups.") end
    this.stopping_powers = stopping_powers
end

"""
    set_momentum_transfer(this::Multigroup_Cross_Sections,momentum_transfer::Vector{Float64})

To set the momentum transfers.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `momentum_transfer::Vector{Float64}` : momentum transfers.

# Output Argument(s)
N/A

"""
function set_momentum_transfer(this::Multigroup_Cross_Sections,momentum_transfer::Vector{Float64})
    if length(momentum_transfer) != this.number_of_groups error("The length of the momentum_transfer don't fit the number of groups.") end
    this.momentum_transfer = momentum_transfer
end

"""
    set_energy_deposition(this::Multigroup_Cross_Sections,
    energy_deposition::Vector{Float64})

To set the energy deposition cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `energy_deposition::Vector{Float64}` : energy deposition cross-sections.

# Output Argument(s)
N/A

"""
function set_energy_deposition(this::Multigroup_Cross_Sections,energy_deposition::Vector{Float64})
    if length(energy_deposition) != this.number_of_groups + 1 error("The length of the energy_deposition cross-sections don't fit the number of groups.") end
    this.energy_deposition = energy_deposition
end

"""
    set_charge_deposition(this::Multigroup_Cross_Sections,
    charge_deposition::Vector{Float64})

To set the charge deposition cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `charge_deposition::Vector{Float64}` : charge deposition cross-sections.

# Output Argument(s)
N/A

"""
function set_charge_deposition(this::Multigroup_Cross_Sections,charge_deposition::Vector{Float64})
    if length(charge_deposition) != this.number_of_groups + 1 error("The length of the charge_deposition cross-sections don't fit the number of groups.") end
    this.charge_deposition = charge_deposition
end

"""
    set_scattering(this::Multigroup_Cross_Sections,scattering::Array{Float64,3})

To set the scattering cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `scattering::Array{Float64,3}` : scattering cross-sections.

# Output Argument(s)
N/A

"""
function set_scattering(this::Multigroup_Cross_Sections,scattering::Array{Float64,3})
    push!(this.scattering,scattering)
end

"""
    get_total(this::Multigroup_Cross_Sections)

Get the total cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `total::Vector{Float64}` : total cross-sections.

"""
function get_total(this::Multigroup_Cross_Sections)
    if ismissing(this.total) error("Unable to get multigroup total cross-sections. Missing data.") end
    return this.total
end

"""
    get_absorption(this::Multigroup_Cross_Sections)

Get the absorption cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `absorption::Vector{Float64}` : absorption cross-sections.

"""
function get_absorption(this::Multigroup_Cross_Sections)
    if ismissing(this.absorption) error("Unable to get multigroup absorption cross-sections. Missing data.") end
    return this.absorption
end

"""
    get_scattering(this::Multigroup_Cross_Sections,index_particle_out::Int64)

Get the scattering cross-sections, for a given outgoing particle.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.
- `index_particle_out::Int64` : index of the outgoing particle.

# Output Argument(s)
- `scattering::Array{Float64}` : scattering cross-sections for a given outgoing particle.

"""
function get_scattering(this::Multigroup_Cross_Sections,index_particle_out::Int64)
    return this.scattering[index_particle_out]
end

"""
    get_boundary_stopping_powers(this::Multigroup_Cross_Sections)

Get the stopping powers at group boundaries.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `boundary_stopping_powers::Vector{Float64}` : stopping powers at group boundaries.

"""
function get_boundary_stopping_powers(this::Multigroup_Cross_Sections)
    return this.boundary_stopping_powers
end

"""
    get_stopping_powers(this::Multigroup_Cross_Sections)

Get the stopping powers.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `stopping_powers::Vector{Float64}` : stopping powers.

"""
function get_stopping_powers(this::Multigroup_Cross_Sections)
    return this.stopping_powers
end

"""
    get_momentum_transfer(this::Multigroup_Cross_Sections)

Get the momentum transfers.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `momentum_transfer::Vector{Float64}` : momentum transfers.

"""
function get_momentum_transfer(this::Multigroup_Cross_Sections)
    return this.momentum_transfer
end

"""
    get_energy_deposition(this::Multigroup_Cross_Sections)

Get the energy deposition cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `energy_deposition::Vector{Float64}` : energy deposition cross-sections.

"""
function get_energy_deposition(this::Multigroup_Cross_Sections)
    return this.energy_deposition
end

"""
    get_charge_deposition(this::Multigroup_Cross_Sections)

Get the charge deposition cross-sections.

# Input Argument(s)
- `this::Multigroup_Cross_Sections` : multigroup cross-sections structure.

# Output Argument(s)
- `charge_deposition::Vector{Float64}` : charge deposition cross-sections.

"""
function get_charge_deposition(this::Multigroup_Cross_Sections)
    return this.charge_deposition
end