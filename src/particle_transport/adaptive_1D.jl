"""
    adaptive_1D(𝒪::Int64,ω::Vector{Float64},𝚽n::Vector{Float64},𝚽12::Float64,s::Float64,Λ::Float64)

Compute the weighting parameters for adaptative calculations over a 1D finite-element. 

# Input Argument(s)
- '𝒪::Int64': order of the flux polynomial expansion.
- 'ω::Vector{Float64}': vector containing the (𝒪+1) weighting factors.
- '𝚽n::Vector{Float64}': moments of the angular flux.
- '𝚽12::Float64': incoming flux.
- 's::Float64': sign corresponding to the sweep direction.
- 'Λ::Float64': incoming/outgoing flux factor ratio.

# Output Argument(s)
- 'isFixed::Bool': indicate if adaptive calculation is completed or not.
- 'ω': corrected weighting factors.

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
function adaptive_1D(𝒪::Int64,ω::Vector{Float64},𝚽n::Vector{Float64},𝚽12::Float64,s::Float64,Λ::Float64)

# Initialization
isFixed = false
ϵ = 1e-16

# Adaptive AWD₀
if 𝒪 == 1

    # Extract weighting factors
    P = -ω[1] / Λ

    if abs(P-1) < ϵ

        # Flux variation in the cell
        u0 = (𝚽12-𝚽n[1])/abs(𝚽n[1])

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
            isFixed = true
        else
            isFixed = false
            ω[1] = -P * Λ
            ω[2] = 1 + P * Λ
        end

    else
        isFixed = true
    end

# Adaptive AWD₁
elseif 𝒪 == 2
    P = ω[1] / Λ
    Q = ω[3] - P * Λ

    if abs(P-0) < ϵ && abs(Q-1) < ϵ

        # Flux variation in the cell
        u0 = (𝚽12-𝚽n[1])
        u1 = -s[1]*sqrt(3)*𝚽n[2]
        if 𝚽n[1] != 0
            u0 = u0/abs(𝚽n[1])
            u1 = u1/abs(𝚽n[1])
        else
            u0 = Inf * sign(u0)
            u1 = Inf * sign(u1)
        end
        if isnan(u0) u0 = 0 end
        if isnan(u1) u1 = 0 end
        s0 = sign(u0); s1 = sign(u1)

        # Weighting factors calculations
        P = 3*(2-3*(1+s1)*abs(u1))/(6+3*(1-s0)*abs(u0)+(6+2*s1)*abs(u1))
        P = min(max(P,0),1)
        if s1 > 0
            Q = ((1-P)-P*abs(u0)*(1-s0)/2-abs(u1)*(19/3-4*P)*P)/(abs(u0)*(2*P^2-5*P+3))
        else
            Q = 1000*abs(u1)/(1+u0)
        end
        Q = min(max(Q,10),1/3)

        # Output data
        if abs(Q - 1) < ϵ && abs(P - 0) < ϵ
            isFixed = true
        else
            isFixed = false
            ω[1] = P * Λ
            ω[2] = 1 - P * Λ
            ω[3] = Q + P * Λ * (2-Q)
        end

    else
        isFixed = true
    end

else
    error("Adaptive scheme is not available for ",𝒪," order flux expansion over the finite-element.")
end

# Return weighting parameters
return isFixed, ω

end