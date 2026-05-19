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
- `Wl::Vector{Wl}` : Legendre moments of the Poškus angular distribution.

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
    if ~haskey(cache_radiant[],"Clk") || length(cache_radiant[]["Clk"][:,1]) < L+1 
    Clk = zeros(L+1,div(L,2)+1)
    for l in range(0,L), k in range(0,div(L,2))
        Clk[l+1,k+1] = (-1)^k * exp( sum(log.(1:2*l-2*k)) - sum(log.(1:k)) - sum(log.(1:l-k)) - sum(log.(1:l-2*k)) )
    end
    cache_radiant[]["Clk"] = Clk
    end

    #----
    # Extract data from cache
    #----
    Clk = cache_radiant[]["Clk"]
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
    Wl = zeros(L+1)
    for l in range(0,L)
        for k in range(0,div(l,2))
            Wlk = 0.0
            Wlk += (A+B)*𝒢a[l-2*k+1]
            for i in range(0,2)
                Wlk += αi[i+1] * 𝒢b[l-2*k+i+1]
            end
            Wl[l+1] += Clk[l+1,k+1] * Wlk
        end
        Wl[l+1] *= 3/(4*(2*A+B)) * (1-C^2)/(2^l)
    end

    #----
    # Correction to deal with high-order Legendre moments
    #----
    for l in range(1,L)
        if abs(Wl[1]) < abs(Wl[l+1])
            Wl[l+1:end] .= 0.0
            break
        end
    end
    return Wl
end