"""
    sauter(Ei::Float64,L::Int64)

Gives the Legendre moments of the Sauter angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wl::Vector{Wl}` : Legendre moments of the Sauter angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled 
  photon-electron-positron cross-sections for the Boltzmann Fokker-Planck equation.

"""
function sauter(Ei::Float64,L::Int64)

    # Compute and save in cache, if not already in cache
    if ~haskey(cache_radiant[],"Clk") || length(cache_radiant[]["Clk"][:,1]) < L+1 
        Clk = zeros(L+1,div(L,2)+1)
        for l in range(0,L), k in range(0,div(L,2))
            Clk[l+1,k+1] = (-1)^k * exp( sum(log.(1:2*l-2*k)) - sum(log.(1:k)) - sum(log.(1:l-k)) - sum(log.(1:l-2*k)) )
        end
        cache_radiant[]["Clk"] = Clk
    end

    # Extract data from cache
    Clk = cache_radiant[]["Clk"]

    # Angular distribution
    Wl = zeros(L+1)
    γ = Ei+1
    β = sqrt(Ei*(Ei+2)/(Ei+1)^2)
    Γ = 1/(4/(3*(1-β^2)^2)+γ*(γ-1)*(γ-2)/(2*β^3)*(2*β/(1-β^2)-log((1+β)/(1-β))))
    α = [1,γ*(γ-1)*(γ-2)/2]
    ΔG3 = zeros(L+3,2)
    for i in range(0,L+2), j in range(0,1)
        ΔG3[i+1,j+1] =  𝒢₃(i,j-4,1,-β,0,1,1)-𝒢₃(i,j-4,1,-β,0,1,-1)
    end
    for l in range(0,L)
        for k in range(0,div(l,2))
            Wlk = 0.0
            for i in range(0,1), j in range(0,1)
                Wlk += α[i+1] * (-1)^j * ΔG3[l-2*k+2*j+1,i+1]
            end
            Wl[l+1] += Clk[l+1,k+1] * Wlk
        end
        Wl[l+1] *= Γ/(2^l)
    end

    # Correction to deal with high-order Legendre moments
    for l in range(1,L)
        if abs(Wl[1]) < abs(Wl[l+1])
            Wl[l+1:end] .= 0.0
            break
        end
    end
    return Wl
end