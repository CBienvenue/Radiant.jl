mutable struct Multigroup_Cross_Sections

    # Variable(s)
    number_of_groups               ::Int64
    total                          ::Union{Missing,Vector{Float64}}
    absorption                     ::Union{Missing,Vector{Float64}}
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
        this.stopping_powers = missing
        this.momentum_transfer = missing
        this.energy_deposition = missing
        this.charge_deposition = missing
        this.scattering = Vector{Array{Float64,3}}()

        return this
    end
end

# Method(s)
Base.propertynames(::Multigroup_Cross_Sections) = 
(
    fieldnames(Multigroup_Cross_Sections)...,
    :set_total,
    :set_absorption,
    :set_stopping_powers,
    :set_momentum_transfer,
    :set_energy_deposition,
    :set_charge_deposition,
    :set_scattering,
    :get_total,
    :get_absorption,
    :get_scattering,
    :get_stopping_powers,
    :get_momentum_transfer,
    :get_energy_deposition,
    :get_charge_deposition
)

function set_total(this::Multigroup_Cross_Sections,total::Vector{Float64})
    if length(total) != this.number_of_groups error("The length of the total cross-sections don't fit the number of groups.") end
    this.total = total
end

function set_absorption(this::Multigroup_Cross_Sections,absorption::Vector{Float64})
    if length(absorption) != this.number_of_groups error("The length of the absorption cross-sections don't fit the number of groups.") end
    this.absorption = absorption
end

function set_stopping_powers(this::Multigroup_Cross_Sections,stopping_powers::Vector{Float64})
    if length(stopping_powers) != this.number_of_groups + 1 error("The length of the stopping_powers don't fit the number of groups.") end
    this.stopping_powers = stopping_powers
end

function set_momentum_transfer(this::Multigroup_Cross_Sections,momentum_transfer::Vector{Float64})
    if length(momentum_transfer) != this.number_of_groups error("The length of the momentum_transfer don't fit the number of groups.") end
    this.momentum_transfer = momentum_transfer
end

function set_energy_deposition(this::Multigroup_Cross_Sections,energy_deposition::Vector{Float64})
    if length(energy_deposition) != this.number_of_groups error("The length of the energy_deposition cross-sections don't fit the number of groups.") end
    this.energy_deposition = energy_deposition
end

function set_charge_deposition(this::Multigroup_Cross_Sections,charge_deposition::Vector{Float64})
    if length(charge_deposition) != this.number_of_groups error("The length of the charge_deposition cross-sections don't fit the number of groups.") end
    this.charge_deposition = charge_deposition
end

function set_scattering(this::Multigroup_Cross_Sections,scattering::Array{Float64,3})
    push!(this.scattering,scattering)
end

function get_total(this::Multigroup_Cross_Sections)
    if ismissing(this.total) error("Unable to get multigroup total cross-sections. Missing data.") end
    return this.total
end

function get_absorption(this::Multigroup_Cross_Sections)
    if ismissing(this.absorption) error("Unable to get multigroup absorption cross-sections. Missing data.") end
    return this.absorption
end

function get_scattering(this::Multigroup_Cross_Sections,index_particle_out::Int64)
    return this.scattering[index_particle_out]
end

function get_stopping_powers(this::Multigroup_Cross_Sections)
    return this.stopping_powers
end

function get_momentum_transfer(this::Multigroup_Cross_Sections)
    return this.momentum_transfer
end

function get_energy_deposition(this::Multigroup_Cross_Sections)
    return this.energy_deposition
end

function get_charge_deposition(this::Multigroup_Cross_Sections)
    return this.charge_deposition
end