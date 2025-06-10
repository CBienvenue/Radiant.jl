"""
    group_structure(Ng::Int64,E::Number,Ec::Number,type::String)

Generate the group structure for multigroup calculations.

# Input Argument(s)
- `Ng::Int64`: number of groups.
- `E::Number`: midpoint energy of the highest energy group.
- `Ec::Number`: cutoff energy.
- `type::String`: type of group structure.

# Output Argument(s)
- `Eᵇ::Vector{Float64}`: vector containing the (Ng+1)-boundaries of the group structure.

# Reference(s)
N/A

"""
function energy_group_structure(Ng::Int64,E::Number,Ec::Number,type::String,custom_energy_boundaries::Union{Vector{Real},Nothing}=nothing)

if type == "linear"

    Emax = E + (E-Ec)/(2*Ng-1)
    Eᵇ = reverse(collect(range(Ec,Emax,Ng+1)))

elseif type == "log"

    f(x) = 2*E - x - 10^(log10(x)-(log10(x)-log10(Ec))/Ng)
    dfdx(x) = -(1+(Ng-1)/Ng*(Ec/x)^(1/Ng))
    x⁻ = 2*E
    x⁺ = Ec
    Emax = newton_bisection(f,dfdx,x⁻,x⁺,1e-7)
    Eᵇ = reverse(10 .^collect(range(log10(Ec),log10(Emax),Ng+1)))

elseif type == "custom"

    if isnothing(custom_energy_boundaries) error("Custom energy boundary vector is empty.") end
    if Ng+1 != length(custom_energy_boundaries) error("The length of custom energy boundary vector is not coherent with the number of energy groups.") end
    if ~all(custom_energy_boundaries[i] > custom_energy_boundaries[i+1] for i in 1:length(custom_energy_boundaries)-1) error("The custom energy boundaries should be strictly decreasing.") end
    Eᵇ = custom_energy_boundaries

else
    error("Unknown group structure.")
end

return Eᵇ
end