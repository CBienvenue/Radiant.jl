"""
    scheme_weights(Nm::Vector{Int64},schemes::Vector{String},Ndims::Int64,isCSD::Bool)

Compute the weights of the closure relations. 

# Input Argument(s)
- `Nm::Vector{Int64}`: vector of orders of the flux polynomial expansion.
- `schemes::Vector{String}`: vector of types of the closure relation.
- `Ndims::Int64`: geometry dimension.
- `isCSD::Bool`: boolean indicating if the continuous slowing-down term is used or not.

# Output Argument(s)
- `ω::Vector{Array{Float64}}`: weighting factors of the scheme.
- `𝒞::Vector{Float64}`: constants related to normalized Legendre expansion.
- `is_adaptive::Bool`: booleans for adaptive calculations.
- `𝒲::Array{Float64}` : weighting constants.

# Reference(s)
- Voloschenko (2011) : Some improvements in solving the transport equation by the use of
  the family of weighted nodal scheme.
- Bienvenue and Hébert (2022) : High-order diamond differencing scheme for the Boltzmann
  Fokker–Planck equation in 1D and 2D Cartesian geometries.

"""
function scheme_weights(Nm::Vector{Int64},schemes::Vector{String},Ndims::Int64,isCSD::Bool)

    # Input validation
    for n in range(1,4)
        if Nm[n] <= 0 error("Scheme required at least a 1st-order polynomial expansion.") end
    end

    # Adaptive keyword
    if any(x->x=="AWD",schemes)
        is_adaptive = true
        schemes .= "AWD"
    else
        is_adaptive = false
    end

    # 1D BTE --------------------------
    if Ndims == 1 && ~isCSD

        xscheme = schemes[1]
        Mx = Nm[1]-1
        ωx = zeros(Mx+2)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            ωx[1] = (-1)^(Mx+1)
            for n in range(0,Mx)
                ωx[n+2] = 1 - (-1)^(Mx+1-n)
            end
        elseif xscheme == "DG" || (xscheme == "AWD" && Mx == 1)
            ωx[1] = 0
            for n in range(0,Mx)
                ωx[n+2] = 1
            end
        elseif xscheme == "DG-"
            if Mx < 1 error("DG- scheme require at least a 2nd-order p expansion.") end
            ωx[1] = 0
            for n in range(0,Mx-1)
                ωx[n+2] = 1
            end
            ωx[Mx+2] = Mx/(2*Mx+1)
        elseif xscheme == "DG+"
            if Mx < 1 error("DG+ scheme require at least a 2nd-order p expansion.") end
            ωx[1] = 0
            for n in range(0,Mx-1)
                ωx[n+2] = 1
            end
            ωx[Mx+2] = 1001
        else
            error("Closure relation not implemented yet.")
        end
        ω = Vector{Array{Float64}}(undef,1)
        ω[1] = ωx

    # 2D BTE --------------------------
    elseif (Ndims == 2 && ~isCSD)
        
        xscheme = schemes[1]
        yscheme = schemes[2]
        Mx = Nm[1]-1
        My = Nm[2]-1
        ωx = zeros(Mx+2,My+1,My+1)
        ωy = zeros(My+2,Mx+1,Mx+1)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            for m in range(0,My)
                ωx[1,m+1,m+1] = (-1)^(Mx+1)
                for n in range(0,Mx)
                    ωx[n+2,m+1,m+1] = 1 - (-1)^(Mx+1-n)
                end
            end
        elseif xscheme == "DG" || (xscheme == "AWD" && Mx == 1)
            for m in range(0,My)
                ωx[1,m+1,m+1] = 0
                for n in range(0,Mx)
                    ωx[n+2,m+1,m+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if yscheme == "DD" || (yscheme == "AWD" && My == 0)
            for m in range(0,Mx)
                ωy[1,m+1,m+1] = (-1)^(My+1)
                for n in range(0,My)
                    ωy[n+2,m+1,m+1] = 1 - (-1)^(My+1-n)
                end
            end
        elseif yscheme == "DG" || (yscheme == "AWD" && My == 1)
            for m in range(0,Mx)
                ωy[1,m+1,m+1] = 0
                for n in range(0,My)
                    ωy[n+2,m+1,m+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        ω = Vector{Array{Float64}}(undef,2)
        ω[1] = ωx
        ω[2] = ωy
        
    # 3D BTE --------------------------
    elseif (Ndims == 3 && ~isCSD)
        
        xscheme = schemes[1]
        yscheme = schemes[2]
        zscheme = schemes[3]
        Mx = Nm[1]-1
        My = Nm[2]-1
        Mz = Nm[3]-1
        ωx = zeros(Mx+2,My+1,Mz+1)
        ωy = zeros(My+2,Mx+1,Mz+1)
        ωz = zeros(Mz+2,Mx+1,My+1)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            for m in range(0,My), t in range(0,Mz)
                ωx[1,m+1,t+1] = (-1)^(Mx+1)
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1] = 1 - (-1)^(Mx+1-n)
                end
            end
        elseif xscheme == "DG"
            for m in range(0,My), t in range(0,Mz)
                ωx[1,m+1,t+1] = 0
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if yscheme == "DD" || (yscheme == "AWD" && My == 0)
            for m in range(0,Mx), t in range(0,Mz)
                ωy[1,m+1,t+1] = (-1)^(My+1)
                for n in range(0,My)
                    ωy[n+2,m+1,t+1] = 1 - (-1)^(My+1-n)
                end
            end
        elseif yscheme == "DG"
            for m in range(0,Mx), t in range(0,Mz)
                ωy[1,m+1,t+1] = 0
                for n in range(0,My)
                    ωy[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if zscheme == "DD" || (zscheme == "AWD" && Mz == 0)
            for m in range(0,Mx), t in range(0,My)
                ωz[1,m+1,t+1] = (-1)^(Mz+1)
                for n in range(0,Mz)
                    ωz[n+2,m+1,t+1] = 1 - (-1)^(Mz+1-n)
                end
            end
        elseif zscheme == "DG"
            for m in range(0,Mx), t in range(0,My)
                ωz[1,m+1,t+1] = 0
                for n in range(0,Mz)
                    ωz[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end

        ω = Vector{Array{Float64}}(undef,3)
        ω[1] = ωx
        ω[2] = ωy
        ω[3] = ωz

    # 1D BFP --------------------------
    elseif (Ndims == 1 && isCSD)
        
        xscheme = schemes[1]
        Escheme = schemes[4]
        Mx = Nm[1]-1
        ME = Nm[4]-1
        ωx = zeros(Mx+2,ME+1,ME+1)
        ωE = zeros(ME+2,Mx+1,Mx+1)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            for m in range(0,ME)
                ωx[1,m+1,m+1] = (-1)^(Mx+1)
                for n in range(0,Mx)
                    ωx[n+2,m+1,m+1] = 1 - (-1)^(Mx+1-n)
                end
            end
        elseif xscheme == "DG" || (xscheme == "AWD" && Mx == 1)
            for m in range(0,ME)
                ωx[1,m+1,m+1] = 0
                for n in range(0,Mx)
                    ωx[n+2,m+1,m+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if Escheme == "DD" || (Escheme == "AWD" && ME == 0)
            for m in range(0,Mx)
                ωE[1,m+1,m+1] = (-1)^(ME+1)
                for n in range(0,ME)
                    ωE[n+2,m+1,m+1] = 1 - (-1)^(ME+1-n)
                end
            end
        elseif Escheme == "DG" || (Escheme == "AWD" && ME == 1)
            for m in range(0,Mx)
                ωE[1,m+1,m+1] = 0
                for n in range(0,ME)
                    ωE[n+2,m+1,m+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        ω = Vector{Array{Float64}}(undef,2)
        ω[1] = ωE
        ω[2] = ωx
    
    # 2D BFP --------------------------
    elseif (Ndims == 2 && isCSD)
        
        xscheme = schemes[1]
        yscheme = schemes[2]
        Escheme = schemes[4]
        Mx = Nm[1]-1
        My = Nm[2]-1
        ME = Nm[4]-1
        ωx = zeros(Mx+2,My+1,ME+1)
        ωy = zeros(My+2,Mx+1,ME+1)
        ωE = zeros(ME+2,Mx+1,My+1)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            for m in range(0,My), t in range(0,ME)
                ωx[1,m+1,t+1] = (-1)^(Mx+1)
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1] = 1 - (-1)^(Mx+1-n)
                end
            end
        elseif xscheme == "DG"
            for m in range(0,My), t in range(0,ME)
                ωx[1,m+1,t+1] = 0
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if yscheme == "DD" || (yscheme == "AWD" && My == 0)
            for m in range(0,Mx), t in range(0,ME)
                ωy[1,m+1,t+1] = (-1)^(My+1)
                for n in range(0,My)
                    ωy[n+2,m+1,t+1] = 1 - (-1)^(My+1-n)
                end
            end
        elseif yscheme == "DG"
            for m in range(0,Mx), t in range(0,ME)
                ωy[1,m+1,t+1] = 0
                for n in range(0,My)
                    ωy[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if Escheme == "DD" || (Escheme == "AWD" && ME == 0)
            for m in range(0,Mx), t in range(0,My)
                ωE[1,m+1,t+1] = (-1)^(ME+1)
                for n in range(0,ME)
                    ωE[n+2,m+1,t+1] = 1 - (-1)^(ME+1-n)
                end
            end
        elseif Escheme == "DG"
            for m in range(0,Mx), t in range(0,My)
                ωE[1,m+1,t+1] = 0
                for n in range(0,ME)
                    ωE[n+2,m+1,t+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end

        ω = Vector{Array{Float64}}(undef,3)
        ω[1] = ωE
        ω[2] = ωx
        ω[3] = ωy

    # 3D BFP --------------------------
    elseif (Ndims == 3 && isCSD)
    
        xscheme = schemes[1]
        yscheme = schemes[2]
        zscheme = schemes[3]
        Escheme = schemes[4]
        Mx = Nm[1]-1
        My = Nm[2]-1
        Mz = Nm[3]-1
        ME = Nm[4]-1
        ωx = zeros(Mx+2,My+1,Mz+1,ME+1)
        ωy = zeros(My+2,Mx+1,Mz+1,ME+1)
        ωz = zeros(Mz+2,Mx+1,My+1,ME+1)
        ωE = zeros(ME+2,Mx+1,My+1,Mz+1)

        if xscheme == "DD" || (xscheme == "AWD" && Mx == 0)
            for m in range(0,My), t in range(0,Mz), w in range(0,ME)
                ωx[1,m+1,t+1,w+1] = (-1)^(Mx+1)
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1,w+1] = 1 - (-1)^(Mx+1-n)
                end
            end
        elseif xscheme == "DG"
            for m in range(0,My), t in range(0,Mz), w in range(0,ME)
                ωx[1,m+1,t+1,w+1] = 0
                for n in range(0,Mx)
                    ωx[n+2,m+1,t+1,w+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if yscheme == "DD" || (yscheme == "AWD" && My == 0)
            for m in range(0,Mx), t in range(0,Mz), w in range(0,ME)
                ωy[1,m+1,t+1,w+1] = (-1)^(My+1)
                for n in range(0,My)
                    ωy[n+2,m+1,t+1,w+1] = 1 - (-1)^(My+1-n)
                end
            end
        elseif yscheme == "DG"
            for m in range(0,Mx), t in range(0,Mz), w in range(0,ME)
                ωy[1,m+1,t+1,w+1] = 0
                for n in range(0,My)
                    ωy[n+2,m+1,t+1,w+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if zscheme == "DD" || (zscheme == "AWD" && Mz == 0)
            for m in range(0,Mx), t in range(0,My), w in range(0,ME)
                ωz[1,m+1,t+1,w+1] = (-1)^(Mz+1)
                for n in range(0,Mz)
                    ωz[n+2,m+1,t+1,w+1] = 1 - (-1)^(Mz+1-n)
                end
            end
        elseif zscheme == "DG"
            for m in range(0,Mx), t in range(0,My), w in range(0,ME)
                ωz[1,m+1,t+1,w+1] = 0
                for n in range(0,Mz)
                    ωz[n+2,m+1,t+1,w+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end
        if Escheme == "DD" || (Escheme == "AWD" && ME == 0)
            for m in range(0,Mx), t in range(0,My), w in range(0,Mz)
                ωE[1,m+1,t+1,w+1] = (-1)^(ME+1)
                for n in range(0,ME)
                    ωE[n+2,m+1,t+1,w+1] = 1 - (-1)^(ME+1-n)
                end
            end
        elseif Escheme == "DG"
            for m in range(0,Mx), t in range(0,My), w in range(0,Mz)
                ωE[1,m+1,t+1,w+1] = 0
                for n in range(0,ME)
                    ωE[n+2,m+1,t+1,w+1] = 1
                end
            end
        else
            error("Closure relation not implemented yet.")
        end

        ω = Vector{Array{Float64}}(undef,4)
        ω[1] = ωE
        ω[2] = ωx
        ω[3] = ωy
        ω[4] = ωz

    else
        error("Wrong number of derivatives.")
    end

    # Constant factors
    𝒞 = zeros(maximum(Nm))
    for i in range(1,maximum(Nm)) 𝒞[i] = sqrt(2*i-1) end
    𝒲 = zeros(Nm[4],Nm[4],Nm[4])
    for i in range(1,Nm[4]), j in range(1,Nm[4]), k in range(1,Nm[4]) 𝒲[i,j,k] = 𝒢₆(i-1,j-1,k-1) end

    return ω, 𝒞, is_adaptive, 𝒲
end

"""
    scheme_weights(Nm::Vector{Int64},schemes::Vector{String},Ndims::Int64,isCSD::Bool)

Compute the weights of the closure relations. 

# Input Argument(s)
- `Nm::Vector{Int64}`: vector of orders of the flux polynomial expansion.
- `schemes::Vector{String}`: vector of types of the closure relation.
- `Ndims::Int64`: geometry dimension.
- `isCSD::Bool`: boolean indicating if the continuous slowing-down term is used or not.

# Output Argument(s)
- `ω::Vector{Array{Float64}}`: weighting factors of the scheme.
- `𝒞::Vector{Float64}`: constants related to normalized Legendre expansion.
- `𝒲::Array{Float64}` : weighting constants.

# Reference(s)
- Voloschenko (2011) : Some improvements in solving the transport equation by the use of
  the family of weighted nodal scheme.
- Bienvenue and Hébert (2022) : High-order diamond differencing scheme for the Boltzmann
  Fokker–Planck equation in 1D and 2D Cartesian geometries.

"""
function scheme_weights_gn(Nm::Vector{Int64},schemes::Vector{String},Ndims::Int64,isCSD::Bool)

    # Input validation
    for n in range(1,4)
        if Nm[n] <= 0 error("Scheme required at least a 1st-order polynomial expansion.") end
    end

    ω = Vector{Vector{Float64}}(undef,4)
    for n in range(1,4)
        ωi = zeros(Nm[n]+1)
        if schemes[n] == "DD"
            ωi[1] = (-1)^Nm[n]
            for m in range(1,Nm[n])
                ωi[m+1] = 1 - (-1)^(Nm[n]-m+1)
            end
        elseif schemes[n] == "DG"
            ωi[1] = 0
            for m in range(1,Nm[n])
                ωi[m+1] = 1
            end
        else
            error("Closure relation not implemented yet.")
        end
        ω[n] = ωi
    end

    # Constant factors
    𝒞 = zeros(maximum(Nm))
    for i in range(1,maximum(Nm)) 𝒞[i] = sqrt(2*i-1) end
    𝒲 = zeros(Nm[4],Nm[4],Nm[4])
    for i in range(1,Nm[4]), j in range(1,Nm[4]), k in range(1,Nm[4]) 𝒲[i,j,k] = 𝒢₆(i-1,j-1,k-1) end

    return ω, 𝒞, 𝒲
end