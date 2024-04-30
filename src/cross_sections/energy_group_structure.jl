"""
    group_structure(Ng::Int64,E::Number,Ec::Number,type::String)

Generate the group structure for multigroup calculations.

# Input Argument(s)
- 'Ng::Int64': number of groups.
- 'E::Number': midpoint energy of the highest energy group.
- 'Ec::Number': cutoff energy.
- 'type::String': type of group structure.

# Output Argument(s)
- 'Eᵇ::Vector{Float64}': vector containing the (Ng+1)-boundaries of the group structure.

# Reference(s)
N/A

"""
function energy_group_structure(Ng::Int64,E::Number,Ec::Number,type::String)

if type == "linear"

    Emax = E + (E-Ec)/(2*Ng-1)
    Eᵇ = reverse(collect(range(Ec,Emax,Ng+1)))

elseif type == "log"

    f(x) = 2*E - x - 10^(log10(x)-(log10(x)-log10(Ec))/Ng)
    dfdx(x) = -(1+(Ng-1)/Ng*(Ec/x)^(1/Ng))
    x⁻ = 2*E
    x⁺ = Ec
    Emax = newton_bissection(f,dfdx,x⁻,x⁺,1e-7)
    Eᵇ = reverse(10 .^collect(range(log10(Ec),log10(Emax),Ng+1)))

else
    error("Unknown group structure.")
end

return Eᵇ
end