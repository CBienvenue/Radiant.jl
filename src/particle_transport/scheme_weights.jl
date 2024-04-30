"""
    scheme_weights(𝒪::Vector{Int64},schemes::Vector{String})

Compute the weights of the closure relations. 

# Input Argument(s)
- '𝒪::Vector{Int64}': vector of orders of the flux polynomial expansion.
- 'schemes::Vector{String}': vector of types of the closure relation.

# Output Argument(s)
- 'ω::Vector{Array{Float64}}': weighting factors of the scheme.
- '𝒞::Vector{Float64}': constants related to normalized Legendre expansion.
- 'is_adaptive::Vector{Bool}': booleans for adaptive calculations.

# Reference(s)
- Voloschenko (2011) : Some improvements in solving the transport equation by the use of
  the family of weighted nodal scheme.
- Bienvenue (2022) : High-order diamond differencing scheme for the Boltzmann Fokker–Planck
  equation in 1D and 2D Cartesian geometries.

"""
function scheme_weights(𝒪::Vector{Int64},schemes::Vector{String})

    𝒞 = Vector{Vector{Float64}}(undef,4)
    ω = Vector{Array{Float64}}(undef,4)
    is_adaptive = zeros(Bool,4)
    𝒪i = [[𝒪[2],𝒪[3],𝒪[4]],[𝒪[1],𝒪[3],𝒪[4]],[𝒪[1],𝒪[2],𝒪[4]],[𝒪[1],𝒪[2],𝒪[3]]]

    # Loop over closure relations
    for n in range(1,4)
        if 𝒪[n] <= 0 error("All scheme require at least a 1st-order p expansion.") end

        # Scheme weighting factors
        if schemes[n] == "DD" || (schemes[n] == "AWD" && 𝒪[n] == 1)
            P = 1
            Q = 1
        elseif schemes[n] == "DG" || (schemes[n] == "AWD" && 𝒪[n] == 2)
            P = 0
            Q = 1
        elseif schemes[n] == "DG-"
            if 𝒪[n] < 2 error("DG- scheme require at least a 2nd-order p expansion.") end
            P = 0
            Q = (𝒪[n] - 1)/(2*𝒪[n] - 1)
        elseif schemes[n] == "DG+"
            if 𝒪[n] < 2 error("DG+ scheme require at least a 2nd-order p expansion.") end
            P = 0
            Q = 300
        else
            error("Closure relation not implemented yet.")
        end

        # Saving weighting factors
        ωi = zeros(𝒪[n]+1,𝒪i[n][1],𝒪i[n][2],𝒪i[n][3])
        for i in range(1,𝒪i[n][1]), j in range(1,𝒪i[n][2]), k in range(1,𝒪i[n][3])
            Pi = Int(i == 1 && j == 1 && k == 1) * P
            if mod(𝒪[n],2) == 1
                ωi[1,i,j,k] = -Pi
                ωi[2:2:𝒪[n]+1,i,j,k] .= 1 + Pi
                if 𝒪[n] >= 2 ωi[3:2:𝒪[n]+1,i,j,k] .= 1 - Pi end
            else
                ωi[1,i,j,k] = Pi
                ωi[2:2:𝒪[n]+1,i,j,k] .= 1 - Pi
                if 𝒪[n] >= 2 ωi[3:2:𝒪[n]+1,i,j,k] .= 1 + Pi end
            end
            ωi[end,i,j,k] = (ωi[end,i,j,k] - 1)*(2-Q) + Q
        end
        ω[n] = ωi

        # Constant factors
        𝒞i = zeros(𝒪[n])
        for i in range(1,𝒪[n]) 𝒞i[i] = sqrt(2*i-1) end
        𝒞[n] = 𝒞i

        # Adaptive keyword
        if (schemes[n] == "AWD") is_adaptive[n] = true end

    end

    return ω, 𝒞, is_adaptive

end