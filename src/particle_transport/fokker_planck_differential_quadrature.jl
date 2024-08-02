"""
    fokker_planck_differential_quadrature(N::Int64,quadrature_type::String,Ndims::Int64)

Calculate the Fokker-Planck scattering matrix using differential quadrature.

# Input Argument(s)
- 'N::Int64': number of directions.
- 'quadrature_type::String': type of quadrature.
- 'Ndims::Int64': dimension of the geometry.

# Output Argument(s)
- 'ℳ::Array{Float64}': Fokker-Planck scattering matrix.
- 'λ₀::Float64': correction factor for total cross-section.

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
    error("Differential quadrature actually unavailable in multidimensional geometry.")
    #=
    Ω₁ = zeros(N,N); I = zeros(N,N); Ω₂ = zeros(N,N); Ω₃ = zeros(N,N); R = zeros(N,N); R1 = zeros(N,N); R2 = zeros(N,N);
    for n in range(1,N)
        Ω₁[n,n] = μ[n]; I[n,n] = 1; Ω₂[n,n] = η[n]; Ω₃[n,n] = ξ[n]; R[n,n] = 1/(1-μ[n]^2); R1 = -ξ[n]/sqrt(1-ξ[n]^2-μ[n]^2); R2 = -η[n]/sqrt(1-η[n]^2-μ[n]^2)
    end 
    
    #----
    # μ-derivative
    #----

    eps = 0.1 # Le choix de ce paramètre est crucial au bon fonctionnement
    
    μi = zeros(0)
    multiplicity = zeros(Int64,0)
    μ_loc = zeros(Int64,N)
    for n in range(1,N)
        if ~any(x->abs(x-μ[n])<eps,μi)
            push!(μi,μ[n])
            push!(multiplicity,0)
        end
        μ_loc[n] = findfirst(x->abs(x-μ[n])<eps,μi)
        multiplicity[μ_loc[n]] += 1 
    end
    Li = length(μi)

    Π = ones(Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Π[n] *= (μi[n]-μi[m]) end
    end

    # First derivative
    Di = zeros(Li,Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,m] = Π[n]/((μi[n]-μi[m])*Π[m]) end
    end
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,n] -= Di[n,m] end
    end

    Dμ = zeros(N,N)
    for n in range(1,N), m in range(1,N)
        Dμ[n,m] = Di[μ_loc[n],μ_loc[m]]/multiplicity[μ_loc[m]]
    end

    # Compute the scattering matrix
    ℳ = zeros(N,N)
    ℳ += (I-Ω₁^2)*Dμ^2 - 2*Ω₁*Dμ
    
    #----
    # η-derivative
    #----

    ηi = zeros(0)
    multiplicity = zeros(Int64,0)
    η_loc = zeros(Int64,N)
    for n in range(1,N)
        if ~any(x->abs(x-η[n])<eps,ηi)
            push!(ηi,η[n])
            push!(multiplicity,0)
        end
        η_loc[n] = findfirst(x->abs(x-η[n])<eps,ηi)
        multiplicity[η_loc[n]] += 1 
    end
    Li = length(ηi)

    Π = ones(Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Π[n] *= (ηi[n]-ηi[m]) end
    end

    # First derivative
    Di = zeros(Li,Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,m] = Π[n]/((ηi[n]-ηi[m])*Π[m]) end
    end
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,n] -= Di[n,m] end
    end

    Dη = zeros(N,N)
    for n in range(1,N), m in range(1,N)
        Dη[n,m] = Di[η_loc[n],η_loc[m]]/multiplicity[η_loc[m]]
    end

    #----
    # ξ-derivative
    #----

    ξi = zeros(0)
    multiplicity = zeros(Int64,0)
    ξ_loc = zeros(Int64,N)
    for n in range(1,N)
        if ~any(x->abs(x-ξ[n])<eps,ξi)
            push!(ξi,ξ[n])
            push!(multiplicity,0)
        end
        ξ_loc[n] = findfirst(x->abs(x-ξ[n])<eps,ξi)
        multiplicity[ξ_loc[n]] += 1 
    end
    Li = length(ξi)

    Π = ones(Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Π[n] *= (ξi[n]-ξi[m]) end
    end

    # First derivative
    Di = zeros(Li,Li)
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,m] = Π[n]/((ξi[n]-ξi[m])*Π[m]) end
    end
    for n in range(1,Li), m in range(1,Li)
        if (n != m) Di[n,n] -= Di[n,m] end
    end

    Dξ = zeros(N,N)
    for n in range(1,N), m in range(1,N)
        Dξ[n,m] = Di[ξ_loc[n],ξ_loc[m]]/multiplicity[ξ_loc[m]]
    end

    if Ndims == 2
        ℳ += R * ( -Ω₂*Dη + Ω₃^2*Dη^2 )
    elseif Ndims == 3
        ℳ += R * ( -Ω₂*Dη - Ω₃*Dξ + Ω₃^2*Dη^2 + Ω₂^2*Dξ^2 - Ω₂*Ω₃*(R1*Dη^2+R2*Dξ^2) )
    end
    
    #----
    # Diagonal positivity
    #----
    
    λ₀ = 0
    for n in range(1,N)
        if (ℳ[n,n] < -λ₀) λ₀ = -ℳ[n,n] end
    end
    for n in range(1,N)
        ℳ[n,n] += λ₀
    end
    =#
end

return ℳ, λ₀

end