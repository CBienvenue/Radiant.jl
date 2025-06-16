"""
    linear_energy_group_structure(Ng::Int64,E::Number,Ec::Number,type::String)

Generate a linear group structure for multigroup calculations.

# Input Argument(s)
- `Ng::Int64`: number of groups.
- `E::Number`: midpoint energy of the highest energy group.
- `Ec::Number`: cutoff energy.

# Output Argument(s)
- `Eᵇ::Vector{Float64}`: vector containing the (Ng+1)-boundaries of the group structure.

# Reference(s)
N/A

"""
function linear_energy_group_structure(Ng::Int64,E::Number,Ec::Number)

    # Validation of input
    if Ec ≥ E error("Cutoff energy is greater than the midpoint energy of the highest energy group.") end
    if Ng < 1 error("Number of groups is less than 1.") end

    # Compute energy boundaries
    Emax = E + (E-Ec)/(2*Ng-1)
    Eᵇ = reverse(collect(range(Ec,Emax,Ng+1)))
    return Eᵇ
end

"""
    log_energy_group_structure(Ng::Int64,E::Number,Ec::Number,type::String)

Generate a logarithmic group structure for multigroup calculations.

# Input Argument(s)
- `Ng::Int64`: number of groups.
- `E::Number`: midpoint energy of the highest energy group.
- `Ec::Number`: cutoff energy.

# Output Argument(s)
- `Eᵇ::Vector{Float64}`: vector containing the (Ng+1)-boundaries of the group structure.

# Reference(s)
N/A

"""
function log_energy_group_structure(Ng::Int64,E::Number,Ec::Number)

    # Validation of input
    if Ec ≥ E error("Cutoff energy is greater than the midpoint energy of the highest energy group.") end
    if Ng < 1 error("Number of groups is less than 1.") end

    # Compute energy boundaries
    f(x) = 2*E - x - 10^(log10(x)-(log10(x)-log10(Ec))/Ng)
    dfdx(x) = -(1+(Ng-1)/Ng*(Ec/x)^(1/Ng))
    x⁻ = 2*E
    x⁺ = Ec
    Emax = newton_bisection(f,dfdx,x⁻,x⁺,1e-7)
    Eᵇ = reverse(10 .^collect(range(log10(Ec),log10(Emax),Ng+1)))
    return Eᵇ
end