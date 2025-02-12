"""
    fokker_planck_differential_quadrature(N::Int64,quadrature_type::String,Ndims::Int64,Qdims::Int64)

Calculate the Fokker-Planck scattering matrix using differential quadrature.

# Input Argument(s)
- 'N::Int64' : number of directions.
- 'quadrature_type::String' : type of quadrature.
- 'Ndims::Int64' : dimension of the geometry.
- 'Qdims::Int64' : dimension of the quadrature.

# Output Argument(s)
- 'ℳ::Array{Float64}' : Fokker-Planck scattering matrix.
- 'λ₀::Float64' : correction factor for total cross-section.

# Reference(s)
- Warsa (2012) : Moment-Preserving SN Discretizations for the One-Dimensional Fokker-Planck
  Equation.

"""
function fokker_planck_differential_quadrature(N::Int64,quadrature_type::String,Ndims::Int64,Qdims::Int64)

if Qdims == 1 

    # Extract quadrature informations
    μ,_ = quadrature(N,quadrature_type,Ndims)

    H = zeros(N,N); I = zeros(N,N)
    for n in range(1,N)
        H[n,n] = μ[n]; I[n,n] = 1
    end 
    
    Π = ones(N)
    for n in range(1,N), m in range(1,N)
        if (n != m) Π[n] *= (μ[n]-μ[m]) end
    end

    # First derivative
    W¹ = zeros(N,N)
    for n in range(1,N), m in range(1,N)
        if (n != m) W¹[n,m] = Π[n]/((μ[n]-μ[m])*Π[m]) end
    end
    for n in range(1,N), m in range(1,N)
        if (n != m) W¹[n,n] -= W¹[n,m] end
    end
    
    # Compute the scattering matrix
    ℳ = zeros(N,N)
    ℳ = (I-H^2)*W¹^2 - 2*H*W¹

    λ₀ = 0
    for n in range(1,N)
        if (ℳ[n,n] < -λ₀) λ₀ = -ℳ[n,n] end
    end
    for n in range(1,N)
        ℳ[n,n] += λ₀
    end

else
    error("Differential quadrature unavailable in multidimensional geometry.")
end

return ℳ, λ₀

end