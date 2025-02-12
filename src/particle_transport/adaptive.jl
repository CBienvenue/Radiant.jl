"""
    adaptive(𝒪x::Int64,ωx::Vector{Float64},hx::Float64,sx::Real,𝚽x12::Float64,
    Qn::Vector{Float64},Σ::Float64)

Compute the weighting parameters for adaptative calculations over a 1D Cartesian
finite-element. 

# Input Argument(s)
- '𝒪x::Int64' : order of the adaptive scheme.
- 'ωx::Vector{Float64}' : weighting factors.
- 'hx::Float64' : inverse of the optical length.
- 'sx::Real' : sign associated with the sweep direction.
- '𝚽x12::Float64' : incoming flux at the edge of the finite-element.
- 'Qn::Vector{Float64}' : source moments.
- 'Σ::Float64' : total cross-section.

# Output Argument(s)
- 'ωx::Vector{Float64}': updated weighting factors.

# Reference(s)
- Alcouffe (1993) : An adaptive weighted diamond differencing method for three-dimensional
  xyz geometry.
- Voloschenko (1994) : Numerical solution of the time-dependant transport equation with
  pulsed sources.
- Germogenova (1994) : Adaptive positive nodal method for transport equation.
- Voloschenko (2011) : Some improvements in solving the transport equation by the use of 
  the family of weighted nodal schemes.

"""
function adaptive(𝒪x::Int64,ωx::Vector{Float64},hx::Float64,sx::Real,𝚽x12::Float64,Qn::Vector{Float64},Σ::Float64)

# Initialization
ϵ = 1e-16

#----
# Adaptive AWD₀
#----
if 𝒪x == 1

    𝚽_DD = (Qn[1] + 2*hx*𝚽x12)/(2*hx+Σ)
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

    ΔQ = sqrt(3) * sx * (Qn[2]-sqrt(3)*hx*sx*𝚽x12)/(Qn[1]+hx*𝚽x12)
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

"""
    adaptive(𝒪x::Int64,𝒪y::Int64,ωx::Array{Float64},ωy::Array{Float64},hx::Float64,
    hy::Float64,sx::Real,sy::Real,𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},
    Qn::Vector{Float64},Σ::Float64)

Compute the weighting parameters for adaptative calculations over a 2D Cartesian
finite-element. 

# Input Argument(s)
- '𝒪x::Int64' : order of the adaptive scheme along x-axis.
- '𝒪y::Int64' : order of the adaptive scheme along y-axis.
- 'ωx::Array{Float64}' : weighting factors along x-axis.
- 'ωy::Array{Float64}' : weighting factors along y-axis.
- 'hx::Float64' : inverse of the optical length along x-axis.
- 'hy::Float64' : inverse of the optical length along y-axis.
- 'sx::Real' : sign associated with the sweep direction along x-axis.
- 'sy::Real' : sign associated with the sweep direction along y-axis.
- '𝚽x12::Float64' : incoming flux at the edge of the finite-element along x-axis.
- '𝚽y12::Float64' : incoming flux at the edge of the finite-element along y-axis.
- 'Qn::Vector{Float64}' : source moments.
- 'Σ::Float64' : total cross-section.

# Output Argument(s)
- 'ωx::Array{Float64}': updated weighting factors along x-axis.
- 'ωy::Array{Float64}': updated weighting factors along y-axis.

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
function adaptive(𝒪x::Int64,𝒪y::Int64,ωx::Array{Float64},ωy::Array{Float64},hx::Float64,hy::Float64,sx::Real,sy::Real,𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},Qn::Vector{Float64},Σ::Float64)

# Initialization    
ϵ = 1e-16

#----
# Adaptive AWD₀
#----
if 𝒪x == 1 && 𝒪y == 1

    𝚽_DD = (Qn[1] + 2*hx*𝚽x12[1] + 2*hy*𝚽y12[1])/(2*hx+2*hy+Σ)
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
            𝚽 = [ Qn[1] + hx*𝚽x12[1] + hy*𝚽y12[1] , Qn[2] - sqrt(3)*hx*sx*𝚽x12[1] + hy*𝚽y12[2] ]
            ΔQ = sx*sqrt(3)*𝚽[2]
        else
            𝚽 = [ Qn[1] + hx*𝚽x12[1] + hy*𝚽y12[1] , Qn[3] + hx*𝚽x12[2] - sqrt(3)*hy*sy*𝚽y12[1] ]
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

"""
    adaptive(𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,ωx::Array{Float64},ωy::Array{Float64},
    ωz::Array{Float64},hx::Float64,hy::Float64,hz::Float64,sx::Real,sy::Real,sz::Real,
    𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},Qn::Vector{Float64},
    Σ::Float64)

Compute the weighting parameters for adaptative calculations over a 3D Cartesian
finite-element. 

# Input Argument(s)
- '𝒪x::Int64' : order of the adaptive scheme along x-axis.
- '𝒪y::Int64' : order of the adaptive scheme along y-axis.
- '𝒪z::Int64' : order of the adaptive scheme along z-axis.
- 'ωx::Vector{Float64}' : weighting factors along x-axis.
- 'ωy::Vector{Float64}' : weighting factors along y-axis.
- 'ωz::Vector{Float64}' : weighting factors along z-axis.
- 'hx::Float64' : inverse of the optical length along x-axis.
- 'hy::Float64' : inverse of the optical length along y-axis.
- 'hz::Float64' : inverse of the optical length along z-axis.
- 'sx::Real' : sign associated with the sweep direction along x-axis.
- 'sy::Real' : sign associated with the sweep direction along y-axis.
- 'sz::Real' : sign associated with the sweep direction along z-axis.
- '𝚽x12::Float64' : incoming flux at the edge of the finite-element along x-axis.
- '𝚽y12::Float64' : incoming flux at the edge of the finite-element along y-axis.
- '𝚽z12::Float64' : incoming flux at the edge of the finite-element along z-axis.
- 'Qn::Vector{Float64}' : source moments.
- 'Σ::Float64' : total cross-section.

# Output Argument(s)
- 'ωx::Vector{Float64}': updated weighting factors along x-axis.
- 'ωy::Vector{Float64}': updated weighting factors along y-axis.
- 'ωz::Vector{Float64}': updated weighting factors along z-axis.

# Reference(s)
- Alcouffe (1993) : An adaptive weighted diamond differencing method for three-dimensional
  xyz geometry.
- Voloschenko (1994) : Numerical solution of the time-dependant transport equation with
  pulsed sources.

"""
function adaptive(𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},hx::Float64,hy::Float64,hz::Float64,sx::Real,sy::Real,sz::Real,𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},Qn::Vector{Float64},Σ::Float64)

# Initialization    
ϵ = 1e-16

#----
# Adaptive AWD₀
#----
if 𝒪x == 1 && 𝒪y == 1 && 𝒪z == 1

    𝚽_DD = (Qn[1] + 2*hx*𝚽x12[1] + 2*hy*𝚽y12[1] + 2*hz*𝚽z12[1])/(2*hx+2*hy+2*hz+Σ)
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

"""
    adaptive(𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,𝒪w::Int64,ωx::Array{Float64},ωy::Array{Float64},
    ωz::Array{Float64},ωw::Array{Float64},hx::Float64,hy::Float64,hz::Float64,hw::Float64,
    sx::Real,sy::Real,sz::Real,sw::Real,𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},
    𝚽z12::Vector{Float64},𝚽w12::Vector{Float64},Qn::Vector{Float64},Σ::Float64)

Compute the weighting parameters for adaptative calculations over a 4D Cartesian
finite-element. 

# Input Argument(s)
- '𝒪x::Int64' : order of the adaptive scheme along x-axis.
- '𝒪y::Int64' : order of the adaptive scheme along y-axis.
- '𝒪z::Int64' : order of the adaptive scheme along z-axis.
- '𝒪w::Int64' : order of the adaptive scheme along w-axis.
- 'ωx::Vector{Float64}' : weighting factors along x-axis.
- 'ωy::Vector{Float64}' : weighting factors along y-axis.
- 'ωz::Vector{Float64}' : weighting factors along z-axis.
- 'ωw::Vector{Float64}' : weighting factors along w-axis.
- 'hx::Float64' : inverse of the optical length along x-axis.
- 'hy::Float64' : inverse of the optical length along y-axis.
- 'hz::Float64' : inverse of the optical length along z-axis.
- 'hw::Float64' : inverse of the optical length along w-axis.
- 'sx::Real' : sign associated with the sweep direction along x-axis.
- 'sy::Real' : sign associated with the sweep direction along y-axis.
- 'sz::Real' : sign associated with the sweep direction along z-axis.
- 'sw::Real' : sign associated with the sweep direction along w-axis.
- '𝚽x12::Float64' : incoming flux at the edge of the finite-element along x-axis.
- '𝚽y12::Float64' : incoming flux at the edge of the finite-element along y-axis.
- '𝚽z12::Float64' : incoming flux at the edge of the finite-element along z-axis.
- '𝚽w12::Float64' : incoming flux at the edge of the finite-element along w-axis.
- 'Qn::Vector{Float64}' : source moments.
- 'Σ::Float64' : total cross-section.

# Output Argument(s)
- 'ωx::Vector{Float64}': updated weighting factors along x-axis.
- 'ωy::Vector{Float64}': updated weighting factors along y-axis.
- 'ωz::Vector{Float64}': updated weighting factors along z-axis.
- 'ωw::Vector{Float64}': updated weighting factors along w-axis.

# Reference(s)
- Alcouffe (1993) : An adaptive weighted diamond differencing method for three-dimensional
  xyz geometry.
- Voloschenko (1994) : Numerical solution of the time-dependant transport equation with
  pulsed sources.

"""
function adaptive(𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,𝒪w::Int64,ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},ωw::Array{Float64},hx::Float64,hy::Float64,hz::Float64,hw::Float64,sx::Real,sy::Real,sz::Real,sw::Real,𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝚽w12::Vector{Float64},Qn::Vector{Float64},Σ::Float64)

    # Initialization    
    ϵ = 1e-16
    
    #----
    # Adaptive AWD₀
    #----
    if 𝒪x == 1 && 𝒪y == 1 && 𝒪z == 1 && 𝒪w == 1

        𝚽_DD = (Qn[1] + 2*hx*𝚽x12[1] + 2*hy*𝚽y12[1] + 2*hz*𝚽z12[1] + 2*hw*𝚽w12[1])/(2*hx+2*hy+2*hz+2*hw+Σ)
        ux = (𝚽x12[1]-𝚽_DD)/𝚽_DD
        uy = (𝚽y12[1]-𝚽_DD)/𝚽_DD
        uz = (𝚽z12[1]-𝚽_DD)/𝚽_DD
        uw = (𝚽w12[1]-𝚽_DD)/𝚽_DD
    
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
        if b*abs(uw) <= 1
            Pw = 1
        else
            Pw = 1/(b*abs(uw))
        end
        if (Px < 0 || isnan(Px)) Px = 0 end
        if (Py < 0 || isnan(Py)) Py = 0 end
        if (Pz < 0 || isnan(Pz)) Pz = 0 end
        if (Pw < 0 || isnan(Pw)) Pw = 0 end
    
        # For 𝚽x12 - (0)
        ωx[1] = -Px      # 𝚽x12 [0]
        ωx[2] = 1 + Px   # 𝚽    [0]
    
        # For 𝚽y12 - (0)
        ωy[1] = -Py      # 𝚽y12 [0]
        ωy[2] = 1 + Py   # 𝚽    [0]
    
        # For 𝚽z12 - (0)
        ωz[1] = -Pz      # 𝚽z12 [0]
        ωz[2] = 1 + Pz   # 𝚽    [0]

        # For 𝚽w12 - (0)
        ωw[1] = -Pw      # 𝚽w12 [0]
        ωw[2] = 1 + Pw   # 𝚽    [0]
    
    else
        error("Adaptive scheme is not available for (𝒪x,𝒪y,𝒪z,𝒪w)=($(𝒪x),$(𝒪y),$(𝒪z),$(𝒪w)) order flux expansion over the finite-element.")
    end
    
    return ωx,ωy,ωz,ωw
end