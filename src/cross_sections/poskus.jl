"""
    poskus(Z::Int64,Ei::Float64,Eγ::Float64,L::Int64)

Gives the Legendre moments of the Bremsstrahlung angular distribution based on the shape
functions of Poškus.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `Eγ::Float64` : Bremsstrahlung photon energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wℓ::Vector{Wℓ}` : Legendre moments of the Poškus angular distribution.

# Reference(s)
- Poškus (2019), Shape functions and singly differential cross sections of bremsstrahlung
  at electron energies from 10 eV to 3 MeV for Z = 1–100.
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled 
  photon-electron-positron cross-sections for the Boltzmann Fokker-Planck equation.

"""
function poskus(Z::Int64,Ei::Float64,Eγ::Float64,L::Int64)

    #----
    # Compute and save in cache, if not already in cache
    #----
    if ~haskey(cache_radiant[],"Cℓk") || length(cache_radiant[]["Cℓk"][:,1]) < L+1 
    Cℓk = zeros(L+1,div(L,2)+1)
    for ℓ in range(0,L), k in range(0,div(L,2))
        Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
    end
    cache_radiant[]["Cℓk"] = Cℓk
    end

    #----
    # Extract data from cache
    #----
    Cℓk = cache_radiant[]["Cℓk"]
    data = fast_load("bremsstrahlung_photons_distribution_poskus_2019.jld2")
    A = data["A"]; B = data["B"]; C = data["C"]
    E = data["E"]; r = data["r"]

    #----
    # Interpolation of parameters A, B and C
    #----
    mₑc² = 0.510999
    if Ei ≥ 3/mₑc²
        β = sqrt(Ei*(Ei+2))/(Ei+1)
        A = 1.0
        B = 0.0
        C = β
    else
        ri = Eγ/Ei
        i = searchsortedfirst(E[Z],Ei)
        i⁻ = searchsortedfirst(r[Z][i-1,:],ri)
        i⁺ = searchsortedfirst(r[Z][i,:],ri)
        if i⁻ < 13
            A⁻ = linear_interpolation(ri,r[Z][i-1,i⁻-1:i⁻],A[Z][i-1,i⁻-1:i⁻])
            B⁻ = linear_interpolation(ri,r[Z][i-1,i⁻-1:i⁻],B[Z][i-1,i⁻-1:i⁻])
            C⁻ = linear_interpolation(ri,r[Z][i-1,i⁻-1:i⁻],C[Z][i-1,i⁻-1:i⁻])
        else
            A⁻ = A[Z][i-1,end]
            B⁻ = B[Z][i-1,end]
            C⁻ = C[Z][i-1,end]
        end
        if i⁺ < 13
            A⁺ = linear_interpolation(ri,r[Z][i,i⁺-1:i⁺],A[Z][i,i⁺-1:i⁺])
            B⁺ = linear_interpolation(ri,r[Z][i,i⁺-1:i⁺],B[Z][i,i⁺-1:i⁺])
            C⁺ = linear_interpolation(ri,r[Z][i,i⁺-1:i⁺],C[Z][i,i⁺-1:i⁺])
        else
            A⁺ = A[Z][i,end]
            B⁺ = B[Z][i,end]
            C⁺ = C[Z][i,end]
        end
        A = log_interpolation(Ei,E[Z][i-1:i],[A⁻,A⁺])
        B = log_interpolation(Ei,E[Z][i-1:i],[B⁻,B⁺])
        C = log_interpolation(Ei,E[Z][i-1:i],[C⁻,C⁺])
    end

    #----
    # Compute Legendre moments of the Poskus angular distribution
    #----
    αi = [C^2,-2*C,1] .* (A-B)
    𝒢a = zeros(L+3)
    𝒢b = zeros(L+3)
    for i in range(0,L+2)
        𝒢a[i+1] = 𝒢₃(i,-2,1,-C,0,1,1)-𝒢₃(i,-2,1,-C,0,1,-1)
        𝒢b[i+1] = 𝒢₃(i,-4,1,-C,0,1,1)-𝒢₃(i,-4,1,-C,0,1,-1)
    end
    Wℓ = zeros(L+1)
    for ℓ in range(0,L)
        for k in range(0,div(ℓ,2))
            Wℓk = 0.0
            Wℓk += (A+B)*𝒢a[ℓ-2*k+1]
            for i in range(0,2)
                Wℓk += αi[i+1] * 𝒢b[ℓ-2*k+i+1]
            end
            Wℓ[ℓ+1] += Cℓk[ℓ+1,k+1] * Wℓk
        end
        Wℓ[ℓ+1] *= 3/(4*(2*A+B)) * (1-C^2)/(2^ℓ)
    end

    #----
    # Correction to deal with high-order Legendre moments
    #----
    for ℓ in range(1,L)
        if abs(Wℓ[1]) < abs(Wℓ[ℓ+1])
            Wℓ[ℓ+1:end] .= 0.0
            break
        end
    end
    return Wℓ
end