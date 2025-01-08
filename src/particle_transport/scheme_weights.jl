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
function scheme_weights(𝒪::Vector{Int64},schemes::Vector{String},Ndims::Int64,isCSD::Bool)

    # Input validation
    for n in range(1,4)
        if 𝒪[n] <= 0 error("All scheme require at least a 1st-order p expansion.") end
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
        Mx = 𝒪[1]-1
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
        ω = ωx

    # 2D BTE --------------------------
    elseif (Ndims == 2 && ~isCSD)
        
        xscheme = schemes[1]
        yscheme = schemes[2]
        Mx = 𝒪[1]-1
        My = 𝒪[2]-1
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
        Mx = 𝒪[1]-1
        My = 𝒪[2]-1
        Mz = 𝒪[3]-1
        ωx = zeros(Mx+2,My+1,Mz+1)
        ωy = zeros(My+2,Mx+1,Mz+1)
        ωz = zeros(My+2,Mx+1,My+1)

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
        Mx = 𝒪[1]-1
        ME = 𝒪[4]-1
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
        error()

    # 3D BFP --------------------------
    elseif (Ndims == 3 && isCSD)
        error()
    else
        error("Wrong number of derivatives.")
    end

    # Constant factors
    𝒞 = zeros(maximum(𝒪))
    for i in range(1,maximum(𝒪)) 𝒞[i] = sqrt(2*i-1) end

    return ω, 𝒞, is_adaptive
end