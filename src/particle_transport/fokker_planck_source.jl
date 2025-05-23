"""
    FP_source(N::Int64,P::Int64,Nm::Int64,T::Vector{Float64},𝚽ℓ::Array{Float64},
    Qℓ::Array{Float64},Ns::Vector{Int64},mat::Array{Int64,3},ℳ::Array{Float64,2},
    Mn::Array{Float64,2},Dn::Array{Float64,2})

Calculate the angular Fokker-Planck source term in Cartesian geometry.

# Input Argument(s)
- `P::Int64`: number of angular interpolation basis.
- `Nm::Int64`: total number of spatial and energy moments.
- `T::Vector{Float64}`: restricted momentum transfer.
- `𝚽ℓ::Array{Float64}`: Legendre components of the in-cell flux.
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.
- `Ns::Vector{Int64}`: number of voxels per axis.    
- `mat::Array{Int64,3}`: material identifier per voxel.
- `ℳ::Array{Float64,2}`: Fokker-Planck scattering matrix.

# Output Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
- Morel (1988) : A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.

"""
function fokker_planck_source(P::Int64,Nm::Int64,T::Vector{Float64},𝚽ℓ::Array{Float64},Qℓ::Array{Float64},Ns::Vector{Int64},mat::Array{Int64,3},ℳ::Array{Float64,2})

# Compute the angular Fokker-Planck source term
for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
  for is in range(1,Nm), n in range(1,P), m in range(1,P)
    Qℓ[n,is,ix,iy,iz] += T[mat[ix,iy,iz]] * ℳ[n,m] * 𝚽ℓ[m,is,ix,iy,iz]
  end
end

return Qℓ
end