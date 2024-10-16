"""
    scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},ndims::Int64,
    Σs::Array{Float64},mat::Array{Int64},type::Int64,P::Int64,pℓ::Vector{Int64},Nm::Int64,
    Ns::Vector{Int64},gi::UnitRange{Int64}=0:0,gf::UnitRange{Int64}=0:0)

Compute the scattering source terms of the transport equation.

# Input Argument(s)
- 'Qℓ::Array{Float64}': Legendre components of the in-cell source.
- '𝚽ℓ::Array{Float64}': Legendre components of the in-cell flux.
- 'ndims::Int64': dimension of the geometry.
- 'Σs::Array{Float64}': Legendre moments of the scattering differential cross-sections.
- 'mat::Array{Int64}': material identifier per voxel.
- 'type::Int64': type of scattering, (1) in-group scattering, (2) out-of-group, scattering,
   (3) from one particle group to another particle one.
- 'P::Int64': number of angular interpolation basis.
- 'pℓ::Vector{Int64}': legendre order associated with each interpolation basis. 
- 'Nm::Int64': number of spatial and/or energy moments.
- 'Ns::Vector{Int64}': number of voxels per axis.
- 'gi::UnitRange{Int64}=0:0': indexes of the initial groups.
- 'gf::UnitRange{Int64}=0:0': indexes of the final groups.

# Output Argument(s)
- 'Qℓ::Array{Float64}': Legendre components of the in-cell source.

# Reference(s)
N/A

"""

# In-scattering
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64})
    @inbounds for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
        Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],pℓ[p]+1] * 𝚽ℓ[p,is,ix,iy,iz]
    end
    return Qℓ 
end

# Out-scattering
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,gf::Int64)
    @inbounds for gi in range(1,Ngi)
        if gi != gf
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
                Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,pℓ[p]+1] * 𝚽ℓ[gi,p,is,ix,iy,iz]
            end
        end
    end
    return Qℓ 
end

# From secondary particle
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,Ngf::Int64,𝚽cutoff::Array{Float64})
    @inbounds for gf in range(1,Ngf), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
        for gi in range(1,Ngi), is in range(1,Nm), p in range(1,P)
            Qℓ[gf,p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,gf,pℓ[p]+1] * 𝚽ℓ[gi,p,is,ix,iy,iz]
        end
        Qℓ[gf,1,1,ix,iy,iz] += Σs[mat[ix,iy,iz],Ngi+1,gf,1] * 𝚽cutoff[1,1,ix,iy,iz]
    end
    return Qℓ 
end