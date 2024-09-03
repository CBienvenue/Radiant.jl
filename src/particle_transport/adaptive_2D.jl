"""
    adaptive_2D(𝒪::Vector{Int64},ω::Vector{Matrix{Float64}},𝚽n::Vector{Float64},𝚽12::Vector{Vector{Float64}},s::Vector{Float64},Λ::Vector{Float64},Ti::Vector{Float64})

Compute the weighting parameters for adaptative calculations over a 2D finite-element. 

# Input Argument(s)
- '𝒪::Vector{Int64}': orders of the flux polynomial expansion.
- 'ω::Vector{Vector{Float64}}': vectors containing the (𝒪+1) weighting factors.
- '𝚽n::Vector{Float64}': moments of the angular flux.
- '𝚽12::Vector{Vector{Float64}}': incoming fluxes.
- 's::Vector{Float64}': signs corresponding to the sweep direction.
- 'Λ::Vector{Float64}': incoming/outgoing flux factors ratio.
- 'Ti::Vector{Float64}': high-order weighting factors.

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
function adaptive_2D(𝒪::Vector{Int64},ω::Vector{Matrix{Float64}},𝚽n::Vector{Float64},𝚽12::Vector{Vector{Float64}},s::Vector{Float64},Λ::Vector{Float64},Ti::Vector{Float64},Qn,h,QQ)

# Initialization
isFixed = zeros(Bool,2)
T = zeros(2)
ϵ = 1e-16

# Adaptive AWD₀
if 𝒪[1] == 1 && 𝒪[2] == 1

    # Loop over all dimension
    for i in range(1,2)
        
        P = -ω[i][1,1] / Λ[i]
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
                ω[i][1,1] = -P * Λ[i]
                ω[i][2,1] = 1 + P * Λ[i]
            end
    
        else
            isFixed[i] = true
        end
    end

# Adaptive mixed AWD₀-AWD₁
elseif (𝒪[1] == 1 && 𝒪[2] == 2) || (𝒪[1] == 2 && 𝒪[2] == 1)

    # Loop over all dimension
    for i in range(1,2)
        
        if 𝒪[i] == 1
            P = -ω[i][1,1] / Λ[i]
            if abs(P-1) < ϵ

                # Flux variation in the cell
                u0 = (𝚽12[i][1]-𝚽n[1])/abs(𝚽n[1])

                # Weighting parameter calculation
                b = 3
                if b*abs(u0) <= 1
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
                    ω[i][1,1] = -P * Λ[i]
                    ω[i][2,1] = 1 + P * Λ[i]
                end
        
            else
                isFixed[i] = true
            end

        elseif 𝒪[i] == 2

            P = ω[i][1,1] / Λ[i]
            Q = ω[i][3,1] - P * Λ[i]

            if abs(P-0) < ϵ && abs(Q-1) < ϵ

                # Flux variation in the cell
                u0 = (𝚽12[i][1]-𝚽n[1])
                u1 = -s[i]*sqrt(3)*𝚽n[2]
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
                P = 0
                if s1 > 0
                    Q = min(max(((1-P)-P*abs(u0)*(1-s0)/2-abs(u1)*(19/3-4*P)*P)/(abs(u0)*(2*P^2-5*P+3)),1/3),1)
                else
                    Q = max(min((abs(u0)+100*abs(u1))/3,100),1) 
                end

                # Output data
                if abs(Q - 1) < ϵ && abs(P - 0) < ϵ
                    isFixed[i] = true
                else
                    isFixed[i] = false
                    ω[i][1,1] = P * Λ[i]
                    ω[i][2,1] = 1 - P * Λ[i]
                    ω[i][3,1] = Q + P * Λ[i] * (2-Q)
                end

            else
                isFixed[i] = true
            end

        end
    end

# Adaptive AWD₁
elseif 𝒪[1] == 2 && 𝒪[2] == 2

    𝚽 = zeros(2)
    x_CM = zeros(2)

    # Loop over all dimension
    for i in range(1,2)

        P = ω[i][1,1] / Λ[i]
        Q = ω[i][3,1] - P * Λ[i]
        T = Ti[i]
        if abs(P-0) < ϵ && abs(Q-1) < ϵ && abs(T-0) < ϵ

            #if (i == 1) 𝚽 = 𝚽n[1:1:𝒪[1]] elseif (i == 2) 𝚽 = 𝚽n[1:𝒪[1]:𝒪[1]*𝒪[2]] end

            if i == 1
                𝚽 = [ Qn[1] + h[1]*𝚽12[1][1] + h[2]*𝚽12[2][1] , Qn[2] - sqrt(3)*h[1]*s[1]*𝚽12[1][1] + h[2]*𝚽12[2][2] ]
            else
                𝚽 = [ Qn[1] + h[1]*𝚽12[1][1] + h[2]*𝚽12[2][1] , Qn[3] + h[1]*𝚽12[1][2] - sqrt(3)*h[2]*s[2]*𝚽12[2][1] ]
            end

            # Flux variation in the cell
            u1 = s[i]*sqrt(3)*𝚽[2]
            if 𝚽[1] != 0
                u1 = u1/𝚽[1]
            elseif abs(u1) < 1e-8
                u1 = 0
            else
                u1 = Inf * sign(u1)
            end
            if isnan(u1) u1 = 0 end

            if u1 > 2.999 u1 = 2.999 end
            if u1 < -2.999 u1 = -2.999 end

            # Estimate the centroid positions
            x_CM_max = 49/100
            u1_temp = copy(u1)
            u1 += 3/2*sign(u1_temp)
            if u1 > 3/2
                if u1 > 3
                    x_CM[i] = x_CM_max
                else
                    x_CM[i] = min(x_CM_max, (2*u1-3+sqrt(12*u1^2-27))/(4*(3*u1)) )
                end
            end
            u1 = u1_temp
            u1 += sign(u1_temp)
            if u1 < -1
                x_CM[i] = max(-x_CM_max,(1+u1)/4)
            end

            if x_CM[i] < -x_CM_max  x_CM[i] = -x_CM_max end
            if x_CM[i] > x_CM_max x_CM[i] = x_CM_max end

            #=
            u1 += sign(u1_temp)
            if u1 > 1
                if u1 < 3
                    Qm = max(1/(u1^u1),1/3)
                else
                    Qm = 1/3
                end
                x_CM[i] = max((1-1/Qm)/4,-x_CM_max)
            end
            u1 += 0.5*sign(u1_temp)
            if u1 < -1
                if u1 > -3
                    Qm = min((-u1)^(-u1),9)
                else
                    Qm = 9
                end
                if Qm < 4/3
                    x_CM[i] = min((Qm-1)/2,x_CM_max)
                else
                    x_CM[i] = min((3*Qm-2+sqrt(9*Qm^2-12*Qm+1))/(18*Qm),x_CM_max)
                end
            end
            =#
        else
            isFixed[i] = true
        end
    end

    # Weighting parameter calculation
    T = zeros(2)
    if ~all(isFixed) 
        P = zeros(2)
        Q,T = constant_linear(x_CM[1],x_CM[2])
        #T *= prod(s)

        # Output data
        for i in range(1,2)
            if abs(Q[i] - 1) < ϵ && abs(P[i] - 0) < ϵ && abs(T[i]) < ϵ
                isFixed[i] = true
            else
                isFixed[i] = false
                ω[i][1,1] = P[i] * Λ[i]
                ω[i][2,1] = 1 - P[i] * Λ[i]
                ω[i][3,1] = Q[i] + P[i] * Λ[i] * (2-Q[i])
                ω[i][1,2] = 0
                ω[i][2,2] = 1
                ω[i][3,2] = Q[i]
            end
        end
    end

else
    error("Adaptive scheme is not available for ",𝒪," order flux expansion over the finite-element.")
end

return all(isFixed), ω, T

end