"""
    scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64})

Compute the elastic (in-group) scattering source.

# Input Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽ℓ::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pℓ::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.

# Output Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64})
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
        Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],pℓ[p]+1] * 𝚽ℓ[p,is,ix,iy,iz]
    end
    return Qℓ 
end

"""
    scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,
    gf::Int64)

Compute the inelastic (out-of-group) scattering source.

# Input Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽ℓ::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pℓ::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Ngi::Int64`: number of energy groups.
- `gf::Int64`: group in which the particles scatter.

# Output Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function scattering_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,gf::Int64)
    for gi in range(1,Ngi)
        if gi != gf
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
                Qℓ[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,pℓ[p]+1] * 𝚽ℓ[gi,p,is,ix,iy,iz]
            end
        end
    end
    return Qℓ 
end

"""
    particle_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,
    Ngf::Int64)

Compute the source produced by a secondary particle.

# Input Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽ℓ::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pℓ::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Ngi::Int64`: number of energy groups for the incoming particle.
- `Ngf::Int64`: number of energy groups for the outgoing particle.

# Output Argument(s)
- `Qℓ::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function particle_source(Qℓ::Array{Float64},𝚽ℓ::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pℓ::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,Ngf::Int64)
    for gf in range(1,Ngf), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
        for gi in range(1,Ngi), is in range(1,Nm), p in range(1,P)
            Qℓ[gf,p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,gf,pℓ[p]+1] * 𝚽ℓ[gi,p,is,ix,iy,iz]
        end
    end
    return Qℓ 
end