"""
    adaptive_3D(𝒪::Vector{Int64},ω::Vector{Matrix{Float64}},𝚽n::Vector{Float64},𝚽12::Vector{Vector{Float64}},Λ::Vector{Float64})

Compute the weighting parameters for adaptative calculations over a 2D finite-element. 

# Input Argument(s)
- '𝒪::Vector{Int64}': orders of the flux polynomial expansion.
- 'ω::Vector{Vector{Float64}}': vectors containing the (𝒪+1) weighting factors.
- '𝚽n::Vector{Float64}': moments of the angular flux.
- '𝚽12::Vector{Vector{Float64}}': incoming fluxes.
- 'Λ::Vector{Float64}': incoming/outgoing flux factors ratio.

# Output Argument(s)
- 'isFixed::Bool': indicate if adaptive calculation is completed or not.
- 'ω': corrected weighting factors.
- 'T': corrected weighting factors.

# Reference(s)
- Alcouffe (1993) : An adaptive weighted diamond differencing method for three-dimensional
  xyz geometry.
- Voloschenko (1994) : Numerical solution of the time-dependant transport equation with
  pulsed sources.
- Germogenova (1994) : Adaptive positive nodal method for transport equation.
- Voloschenko (2011) : Some improvements in solving the transport equation by the use of 
  the family of weighted nodal schemes.
- Bienvenue (2023) : Adaptive Gradient-Driven Coupled Linear Schemes and their Usefulness
  for Charged Particle Transport.

"""
function adaptive_3D(𝒪::Vector{Int64},ω,𝚽n::Vector{Float64},𝚽12,Λ::Vector{Float64})

# Initialization
isFixed = zeros(Bool,3)
ϵ = 1e-16

# Adaptive AWD₀
if 𝒪[1] == 1 && 𝒪[2] == 1 && 𝒪[3] == 1

    # Loop over all dimension
    for i in range(1,3)
        
        P = -ω[i][1,1,1] / Λ[i]
        if abs(P-1) < ϵ

            # Flux variation in the cell
            u0 = (𝚽12[i][1]-𝚽n[1])/abs(𝚽n[1])
    
            # Weighting parameter calculation
            b = 3
            if u0 < 0 || b*abs(u0) <= 1
                P = 1
            else
                P = 1/(b*abs(u0))
            end
            if (P < 0 || isnan(P)) P = 0 end
            
            # Output data
            if abs(P-1) < ϵ
                isFixed[i] = true
            else
                isFixed[i] = false
                ω[i][1,1,1] = -P * Λ[i]
                ω[i][2,1,1] = 1 + P * Λ[i]
            end
    
        else
            isFixed[i] = true
        end
    end

else
    error("Adaptive scheme is not available for ",𝒪," order flux expansion over the finite-element.")
end

return all(isFixed), ω

end