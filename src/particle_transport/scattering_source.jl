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
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},ndims::Int64,Σs::Array{Float64},mat::Array{Int64},type::Int64,P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},gi::UnitRange{Int64}=0:0,gf::UnitRange{Int64}=0:0)

if type == 1 # In-scattering
    @inbounds for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
        Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],pℓ[p]+1] * 𝚽ℓ[p,is,ix,iy,iz]
    end
elseif type == 2 # Out-scattering
    igf = gf[1]
    @inbounds for igi in gi
        if igi != igf
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
                Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],igi,pℓ[p]+1] * 𝚽ℓ[igi,p,is,ix,iy,iz]
            end
        end
    end
elseif type == 3 # From secondary particle
    @inbounds for igi in gi, igf in gf, ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
        Qℓ[igf,p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],igi,igf,pℓ[p]+1] * 𝚽ℓ[igi,p,is,ix,iy,iz]
    end
else
    error(" Incorrect type of scattering source.")
end

return Qℓ 
end