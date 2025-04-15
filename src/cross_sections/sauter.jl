"""
    sauter(Ei::Float64,L::Int64)

Gives the Legendre moments of the Sauter angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wℓ::Vector{Wℓ}` : Legendre moments of the Sauter angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled 
  photon-electron-positron cross-sections for the Boltzmann Fokker-Planck equation.

"""
function sauter(Ei::Float64,L::Int64)

    # Compute and save in cache, if not already in cache
    if ~haskey(cache_radiant[],"Cℓk") || length(cache_radiant[]["Cℓk"][:,1]) < L+1 
        Cℓk = zeros(L+1,div(L,2)+1)
        for ℓ in range(0,L), k in range(0,div(L,2))
            Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
        end
        cache_radiant[]["Cℓk"] = Cℓk
    end

    # Extract data from cache
    Cℓk = cache_radiant[]["Cℓk"]

    # Angular distribution
    Wℓ = zeros(L+1)
    γ = Ei+1
    β = sqrt(Ei*(Ei+2)/(Ei+1)^2)
    Γ = 1/(4/(3*(1-β^2)^2)+γ*(γ-1)*(γ-2)/(2*β^3)*(2*β/(1-β^2)-log((1+β)/(1-β))))
    α = [1,γ*(γ-1)*(γ-2)/2]
    ΔG3 = zeros(L+3,2)
    for i in range(0,L+2), j in range(0,1)
        ΔG3[i+1,j+1] =  𝒢₃(i,j-4,1,-β,0,1,1)-𝒢₃(i,j-4,1,-β,0,1,-1)
    end
    for ℓ in range(0,L)
        for k in range(0,div(ℓ,2))
            Wℓk = 0.0
            for i in range(0,1), j in range(0,1)
                Wℓk += α[i+1] * (-1)^j * ΔG3[ℓ-2*k+2*j+1,i+1]
            end
            Wℓ[ℓ+1] += Cℓk[ℓ+1,k+1] * Wℓk
        end
        Wℓ[ℓ+1] *= Γ/(2^ℓ)
    end

    # Correction to deal with high-order Legendre moments
    for ℓ in range(1,L)
        if abs(Wℓ[1]) < abs(Wℓ[ℓ+1])
            Wℓ[ℓ+1:end] .= 0.0
            break
        end
    end
    return Wℓ
end