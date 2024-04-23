"""
    FP_scattering_matrix(N::Int64,Ω::Union{Vector{Vector{Float64}},Vector{Float64}},
    w::Vector{Float64},Ndims::Int64,method::String,Mn::Array{Float64},Dn::Array{Float64},
    pℓ::Vector{Int64},P::Int64)

Calculate the Fokker-Planck scattering matrix.

See also [`fokker_planck_source`](@ref), [`angular_polynomial_basis`](@ref).

# Input Argument(s)
- 'N::Int64': number of directions.
- 'Ω::Union{Vector{Vector{Float64}},Vector{Float64}}': director cosines.
- 'w::Vector{Float64}': quadrature weights.
- 'Ndims::Int64': dimension of the geometry.
- 'method::String': discretisation method for the angular Fokker-Planck term.
- 'Mn::Array{Float64}': moment-to-discrete matrix.
- 'Dn::Array{Float64}': discrete-to-moment matrix.
- 'pℓ::Vector{Int64}': legendre order associated with each interpolation basis. 
- 'P::Int64': number of angular interpolation basis.

# Output Argument(s)
- 'ℳ::Array{Float64}': Fokker-Planck scattering matrix.
- 'λ₀::Float64': correction factor for total cross-section.

# Author(s)
Charles Bienvenue

# Reference(s)
- Warsa (2012) : Moment-Preserving SN Discretizations for the One-Dimensional Fokker-Planck
  Equation.

"""
function fokker_planck_scattering_matrix(N::Int64,Nd::Int64,quadrature_type::String,Ndims::Int64,method::String,Mn::Array{Float64},Dn::Array{Float64},pℓ::Vector{Int64},pm::Vector{Int64},P::Int64)

if method == "finite-element"
    ℳ, λ₀, Mn_FP, Dn_FP, N_FP = fokker_planck_finite_element(N,quadrature_type,Ndims,pℓ,pm,P,Nd,Mn,Dn)
elseif method == "differential-quadrature"
    ℳ, λ₀ = fokker_planck_differential_quadrature(N,quadrature_type,Ndims)
    Mn_FP = Mn; Dn_FP = Dn; N_FP = Nd
elseif method == "galerkin"
    ℳ, λ₀ = fokker_planck_galerkin(Nd,Mn,Dn,pℓ,P)
    Mn_FP = Mn; Dn_FP = Dn; N_FP = Nd
else
    error("Unknown method to treat the Fokker-Planck term.")
end

ℳ = Dn_FP * ℳ * Mn_FP

return ℳ, λ₀, Mn_FP, Dn_FP, N_FP
end