"""
    biggs_lighthill(Z::Int64,Ei::Float64,ρ::Float64=density(Z))

Gives the photoelectric absorption cross-section based on the Biggs and Lighthill
interpolation.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.

# Output Argument(s)
- `σa::Float64` : photoelectric differential cross-sections.

# Reference(s)
- Biggs and Lighthill (1988), Analytical Approximations for X-Ray Cross Sections III
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.

"""
function biggs_lighthill(Z::Int64,Ei::Float64)

    # Extract data
    data = fast_load("photoelectric_biggs_lighthill_1988.jld2")

    # Compute interpolation parameters
    E⁻ = data["E"][Z]
    M = data["M"][Z]
    A = Vector{Float64}(undef,4)
    N_interval = length(E⁻)
    for i in range(1,N_interval)
        if i == N_interval && Ei >= E⁻[i] || Ei >= E⁻[i] && Ei < E⁻[i+1] 
            A = M[i,:]
            break
        end
    end

    # Absorption cross-sections
    σa = 0.0
    for i in range(1,4)
        σa += A[i]/Ei^i
    end
    σa *= density(Z) / nuclei_density(Z,density(Z))
    return σa
end