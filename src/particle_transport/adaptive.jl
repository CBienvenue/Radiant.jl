function adaptive(𝒪x::Int64,ωx,μ::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Float64,Σ::Float64)

# Initialization    
ϵ = 1e-16
sx = sign(μ)

#----
# Adaptive AWD₀
#----
if 𝒪x == 1

    hx = abs(μ)/(Σ*Δx)
    𝚽_DD = (Qn[1]/Σ + 2*hx*𝚽x12)/(2*hx+1)
    ux = (𝚽x12-𝚽_DD)/𝚽_DD

    b = 3
    if b*abs(ux) <= 1
        Px = 1
    else
        Px = 1/(b*abs(ux))
    end
    if (Px < 0 || isnan(Px)) Px = 0 end

    # For 𝚽x12 - (0)
    ωx[1] = -Px      # 𝚽x12 [0]
    ωx[2] = 1 + Px   # 𝚽    [0]

#----
# Adaptive AWD₁
#----
elseif 𝒪x == 2

    ΔQ = sqrt(3) * sx * (Qn[2]*Δx-sqrt(3)*abs(μ)*sx*𝚽x12)/(Qn[1]*Δx+abs(μ)*𝚽x12)
    if ΔQ > 3/2
        Qx = (3+ΔQ)/(3*(3-ΔQ))
    elseif ΔQ < -1
        Qx = -1/ΔQ 
    else
        Qx = 1
    end
    Px = 0

    # For 𝚽x12 - (0)
    ωx[1] = Px        # 𝚽x12 [0]
    ωx[2] = 1 - Px    # 𝚽    [0]
    ωx[3] = Qx + Px   # 𝚽    [1]

else
    error("Adaptive scheme is not available for 𝒪x=$(𝒪x) order flux expansion over the finite-element.")
end

return ωx
end

function adaptive(𝒪x,𝒪y,ωx,ωy,μ,η,Δx,Δy,Qn,𝚽x12,𝚽y12,Σ)

# Initialization    
ϵ = 1e-16
sx = sign(μ)
sy = sign(η)

#----
# Adaptive AWD₀
#----
if 𝒪x == 1 && 𝒪y == 1

    hx = abs(μ)/(Σ*Δx)
    hy = abs(η)/(Σ*Δy)
    𝚽_DD = (Qn[1]/Σ + 2*hx*𝚽x12[1] + 2*hy*𝚽y12[1])/(2*hx+2*hy+1)
    ux = (𝚽x12[1]-𝚽_DD)/𝚽_DD
    uy = (𝚽y12[1]-𝚽_DD)/𝚽_DD

    b = 3
    if b*abs(ux) <= 1
        Px = 1
    else
        Px = 1/(b*abs(ux))
    end
    if b*abs(uy) <= 1
        Py = 1
    else
        Py = 1/(b*abs(uy))
    end
    if (Px < 0 || isnan(Px)) Px = 0 end
    if (Py < 0 || isnan(Py)) Py = 0 end

    # For 𝚽x12 - (0)
    ωx[1] = -Px      # 𝚽x12 [0]
    ωx[2] = 1 + Px   # 𝚽    [0]

    # For 𝚽y12 - (0)
    ωy[1] = -Py      # 𝚽y12 [0]
    ωy[2] = 1 + Py   # 𝚽    [0]

#----
# Adaptive AWD₁
#----
elseif 𝒪x == 2 && 𝒪y == 2

    # Loop over both axis to estimate CM position
    x_CM = zeros(2)
    for i in range(1,2)
        if i == 1
            𝚽 = [ Qn[1] + abs(μ)/Δx*𝚽x12[1] + abs(η)/Δy*𝚽y12[1] , Qn[2] - sqrt(3)*abs(μ)/Δx*sx*𝚽x12[1] + abs(η)/Δy*𝚽y12[2] ]
            ΔQ = sx*sqrt(3)*𝚽[2]
        else
            𝚽 = [ Qn[1] + abs(μ)/Δx*𝚽x12[1] + abs(η)/Δy*𝚽y12[1] , Qn[3] + abs(μ)/Δx*𝚽x12[2] - sqrt(3)*abs(η)/Δy*sy*𝚽y12[1] ]
            ΔQ = sy*sqrt(3)*𝚽[2]
        end

        # Flux variation in the cell
        if 𝚽[1] != 0
            ΔQ = ΔQ/𝚽[1]
        elseif abs(ΔQ) < 1e-8
            ΔQ = 0
        else
            ΔQ = Inf * sign(ΔQ)
        end
        if isnan(ΔQ) ΔQ = 0 end
        if ΔQ > 2.999 ΔQ = 2.999 end
        if ΔQ < -2.999 ΔQ = -2.999 end

        # Estimate the centroid positions
        x_CM_max = 49/100
        ΔQ_temp = copy(ΔQ)
        ΔQ += 3/2*sign(ΔQ_temp)
        if ΔQ > 3/2
            if ΔQ > 3
                x_CM[i] = x_CM_max
            else
                x_CM[i] = min(x_CM_max, (2*ΔQ-3+sqrt(12*ΔQ^2-27))/(4*(3*ΔQ)) )
            end
        end
        ΔQ = ΔQ_temp
        ΔQ += sign(ΔQ_temp)
        if ΔQ < -1
            x_CM[i] = max(-x_CM_max,(1+ΔQ)/4)
        end
        if x_CM[i] < -x_CM_max  x_CM[i] = -x_CM_max end
        if x_CM[i] > x_CM_max x_CM[i] = x_CM_max end
    end


    # Compute Qx, Qy, Tx and Ty from CM position
    Q,T = constant_linear(x_CM[1],x_CM[2])

    # For 𝚽x12 - (0)
    ωx[2,1,1] = 1                 # 𝚽    [0,0]
    ωx[3,1,1] = Q[1]              # 𝚽    [1,0]

    # For 𝚽x12 - (1)
    ωx[2,2,2] = 1                 # 𝚽    [0,1]
    ωx[3,2,2] = Q[1]              # 𝚽    [1,1]
    ωx[3,1,2] = sy*T[1]/sqrt(3)   # 𝚽    [1,0]

    # For 𝚽y12 - (0)
    ωy[2,1,1] = 1                 # 𝚽    [0,0]
    ωy[3,1,1] = Q[2]              # 𝚽    [0,1]

    # For 𝚽y12 - (1)
    ωy[2,2,2] = 1                 # 𝚽    [1,0]
    ωy[3,2,2] = Q[2]              # 𝚽    [1,1]
    ωy[3,1,2] = sx*T[2]/sqrt(3)   # 𝚽    [0,1]

else
    error("Adaptive scheme is not available for (𝒪x,𝒪y)=($(𝒪x),$(𝒪y)) order flux expansion over the finite-element.")
end

return ωx,ωy
end

function adaptive(𝒪x,𝒪y,𝒪z,ωx,ωy,ωz,μ,η,ξ,Δx,Δy,Δz,Qn,𝚽x12,𝚽y12,𝚽z12,Σ)

# Initialization    
ϵ = 1e-16
sx = sign(μ)
sy = sign(η)
sz = sign(ξ)

#----
# Adaptive AWD₀
#----
if 𝒪x == 1 && 𝒪y == 1 && 𝒪z == 1

    hx = abs(μ)/(Σ*Δx)
    hy = abs(η)/(Σ*Δy)
    hz = abs(ξ)/(Σ*Δz)
    𝚽_DD = (Qn[1]/Σ + 2*hx*𝚽x12[1] + 2*hy*𝚽y12[1] + 2*hz*𝚽z12[1])/(2*hx+2*hy+2*hz+1)
    ux = (𝚽x12[1]-𝚽_DD)/𝚽_DD
    uy = (𝚽y12[1]-𝚽_DD)/𝚽_DD
    uz = (𝚽z12[1]-𝚽_DD)/𝚽_DD

    b = 3
    if b*abs(ux) <= 1
        Px = 1
    else
        Px = 1/(b*abs(ux))
    end
    if b*abs(uy) <= 1
        Py = 1
    else
        Py = 1/(b*abs(uy))
    end
    if b*abs(uz) <= 1
        Pz = 1
    else
        Pz = 1/(b*abs(uz))
    end
    if (Px < 0 || isnan(Px)) Px = 0 end
    if (Py < 0 || isnan(Py)) Py = 0 end
    if (Pz < 0 || isnan(Pz)) Pz = 0 end

    # For 𝚽x12 - (0)
    ωx[1] = -Px      # 𝚽x12 [0]
    ωx[2] = 1 + Px   # 𝚽    [0]

    # For 𝚽y12 - (0)
    ωy[1] = -Py      # 𝚽y12 [0]
    ωy[2] = 1 + Py   # 𝚽    [0]

    # For 𝚽z12 - (0)
    ωz[1] = -Pz      # 𝚽z12 [0]
    ωz[2] = 1 + Pz   # 𝚽    [0]

else
    error("Adaptive scheme is not available for (𝒪x,𝒪y,𝒪z)=($(𝒪x),$(𝒪y),$(𝒪z)) order flux expansion over the finite-element.")
end

return ωx,ωy,ωz
end