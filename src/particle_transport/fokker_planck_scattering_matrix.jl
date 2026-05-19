"""
    FP_scattering_matrix(N::Int64,Ω::Union{Vector{Vector{Float64}},Vector{Float64}},
    w::Vector{Float64},Ndims::Int64,method::String,Mn::Array{Float64},Dn::Array{Float64},
    pl::Vector{Int64},P::Int64)

Calculate the Fokker-Planck scattering matrix for the discrete ordinates method.

# Input Argument(s)
- `N::Int64`: number of directions.
- `Ω::Union{Vector{Vector{Float64}},Vector{Float64}}`: director cosines.
- `w::Vector{Float64}`: quadrature weights.
- `Ndims::Int64`: dimension of the geometry.
- `method::String`: discretisation method for the angular Fokker-Planck term.
- `Mn::Array{Float64}`: moment-to-discrete matrix.
- `Dn::Array{Float64}`: discrete-to-moment matrix.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `P::Int64`: number of angular interpolation basis.
- `Qdims::Int64`: quadrature dimension.

# Output Argument(s)
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `λ₀::Float64`: correction factor for total cross-section.

# Reference(s)
- Warsa and Prinja (2012) : Moment-Preserving SN Discretizations for the One-Dimensional   
  Fokker-Planck Equation.

"""
function fokker_planck_scattering_matrix(N::Int64,Nd::Int64,quadrature_type::String,Ndims::Int64,method::String,Mn::Array{Float64},Dn::Array{Float64},pl::Vector{Int64},P::Int64,Qdims::Int64)
    if method == "finite-difference"
        ℳ, λ₀ = fokker_planck_finite_difference(N,quadrature_type,Ndims,Qdims,Mn,Dn)
    elseif method == "differential-quadrature"
        ℳ, λ₀ = fokker_planck_differential_quadrature(N,quadrature_type,Ndims,Qdims,Mn,Dn)
    elseif method == "galerkin"
        ℳ, λ₀ = fokker_planck_galerkin(Nd,Mn,Dn,pl,P)
    else
        error("Unknown method to treat the Fokker-Planck term.")
    end
    return ℳ, λ₀
end

"""
    FP_scattering_matrix(method::String,pl::Vector{Int64},P::Int64)

Calculate the Fokker-Planck scattering matrix for the spherical harmonics method.

# Input Argument(s)
- `method::String`: discretisation method for the angular Fokker-Planck term.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `P::Int64`: number of angular interpolation basis.

# Output Argument(s)
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `λ₀::Float64`: correction factor for total cross-section.

# Reference(s)
- Warsa and Prinja (2012) : Moment-Preserving SN Discretizations for the One-Dimensional   
  Fokker-Planck Equation.

"""
function fokker_planck_scattering_matrix(method::String,pl::Vector{Int64},P::Int64)
    if method == "galerkin"
        ℳ = fokker_planck_galerkin(P,pl)
    else
        error("Unknown method to treat the Fokker-Planck term.")
    end
    # Stabilization of the scattering matrix
    λ₀ = 0.0
    for p in range(1,P)
        if (ℳ[p,p] < -λ₀) λ₀ = -ℳ[p,p] end
    end
    for p in range(1,P)
        ℳ[p,p] += λ₀
    end
    return ℳ, λ₀
end