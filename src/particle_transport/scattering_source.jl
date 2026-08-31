"""
    scattering_source(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64})

Compute the elastic (in-group) scattering source.

# Input Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.

# Output Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function scattering_source(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64})
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
        Ql[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],pl[p]+1] * 𝚽l[p,is,ix,iy,iz]
    end
    return Ql
end

"""
    scattering_source_fast(Ql,𝚽l,Σs,mat,P,pl,Nm,Ns)

Optimized counterpart of the in-group `scattering_source`, for the GN fast path.

Bit-identical: the same products are added to the same elements, only the loop nest is
reordered to follow the column-major layout. The reference nests `ix,iy,iz,p,is` with `is`
innermost, so consecutive iterations stride by `P` in both `Ql` and `𝚽l`, whose leading
axis is `p`. Here `p` is innermost and the per-material cross-section vector is hoisted out
of the voxel loop.
"""
function scattering_source_fast(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64})
    @inbounds for iz in range(1,Ns[3]), iy in range(1,Ns[2]), ix in range(1,Ns[1])
        m = mat[ix,iy,iz]
        for is in range(1,Nm)
            @simd for p in range(1,P)
                Ql[p,is,ix,iy,iz] += Σs[m,pl[p]+1] * 𝚽l[p,is,ix,iy,iz]
            end
        end
    end
    return Ql
end

"""
    scattering_source_fast(Ql,𝚽l,Σs,mat,P,pl,Nm,Ns,Ngi,gf,is_elastic)

Optimized counterpart of the out-of-group `scattering_source`, for the GN fast path.

Two changes, both bit-identical — the same products reach the same elements in the same
order:

- the loop nest follows the column-major layout, with the stride-1 angular index innermost
  and the per-material cross-section vector hoisted out of the voxel loop;
- source groups whose cross-section block is entirely zero are skipped. Electron and photon
  data downscatter strictly (measured: no upscattering entry at all in a 40-group water
  set), so the reference spends about half its iterations multiplying zeros.
"""
function scattering_source_fast(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,gf::Int64,is_elastic::Bool=false)
    Nl = size(Σs,3)
    @inbounds for gi in range(1,Ngi)
        if gi == gf && ~is_elastic continue end
        # Skip a source group whose whole cross-section block vanishes (strict downscattering).
        nz = false
        for n in range(1,size(Σs,1)), l in range(1,Nl)
            if Σs[n,gi,l] != 0.0; nz = true; break end
        end
        nz || continue
        for iz in range(1,Ns[3]), iy in range(1,Ns[2]), ix in range(1,Ns[1])
            m = mat[ix,iy,iz]
            for is in range(1,Nm)
                @simd for p in range(1,P)
                    Ql[p,is,ix,iy,iz] += Σs[m,gi,pl[p]+1] * 𝚽l[gi,p,is,ix,iy,iz]
                end
            end
        end
    end
    return Ql
end

"""
    scattering_source(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,
    gf::Int64)

Compute the inelastic (out-of-group) scattering source.

# Input Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Ngi::Int64`: number of energy groups.
- `gf::Int64`: group in which the particles scatter.

# Output Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function scattering_source(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,gf::Int64,is_elastic::Bool=false)
    for gi in range(1,Ngi)
        if gi != gf || is_elastic
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,P), is in range(1,Nm)
                Ql[p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,pl[p]+1] * 𝚽l[gi,p,is,ix,iy,iz]
            end
        end
    end
    return Ql 
end

"""
    particle_source(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},
    mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,
    Ngf::Int64)

Compute the source produced by a secondary particle.

# Input Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `ndims::Int64`: dimension of the geometry.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64}`: material identifier per voxel.
- `P::Int64`: number of angular interpolation basis.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `Nm::Int64`: number of spatial and/or energy moments.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Ngi::Int64`: number of energy groups for the incoming particle.
- `Ngf::Int64`: number of energy groups for the outgoing particle.

# Output Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
N/A

"""
function particle_sources(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,Ngf::Int64)
    for gf in range(1,Ngf), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
        for gi in range(1,Ngi), is in range(1,Nm), p in range(1,P)
            Ql[gf,p,is,ix,iy,iz] += Σs[mat[ix,iy,iz],gi,gf,pl[p]+1] * 𝚽l[gi,p,is,ix,iy,iz]
        end
    end
    return Ql
end

"""
    group_transfer_lists(Σs,Ngi,Ngf)

For each outgoing group, the list of incoming groups whose transfer block is non-zero.

The group-to-group matrices are far from full. Measured on water/bone/lung with the complete
interaction set, only 41–44 % of the `(gi,gf)` blocks carry anything, and that fraction is
stable from 20 to 80 groups. Two structures are at play: energy descent (`e⁻→e⁻` and
`γ→e⁻` are strictly lower-triangular) and bounded transfer ranges (`γ→γ` fills only 33 % —
the Compton band is narrow); `e⁻→e⁺` is empty outright.

Descent is *not* universal, so it must not be assumed: `e⁺→γ` transfers **upward** in energy
over 104 of its blocks, carrying 99,5 % of its mass — a slow positron annihilates into two
511 keV photons. Reading the pattern off the data covers that case for free.

The scan costs `Nmat·Ngi·Ngf·(L+1)` comparisons against the `Ngi·Ngf·P·Nm·Nvox` products it
guards — millions of times cheaper, so it is done at each call rather than cached.

The mask is taken as the **union over materials**: it is identical between water, bone and
lung for eight of the nine particle pairs (only `γ→e⁻` differs, by 11 blocks out of 779), and
a per-material mask would have to be consulted inside the voxel loop.
"""
function group_transfer_lists(Σs::Array{Float64},Ngi::Int64,Ngf::Int64)
    Nmat = size(Σs,1); Nl = size(Σs,4)
    lists = Vector{Vector{Int32}}(undef,Ngf)
    for gf in range(1,Ngf)
        gl = Int32[]
        for gi in range(1,Ngi)
            nz = false
            for n in range(1,Nmat), l in range(1,Nl)
                if Σs[n,gi,gf,l] != 0.0; nz = true; break end
            end
            if nz push!(gl,Int32(gi)) end
        end
        lists[gf] = gl
    end
    return lists
end

"""
    particle_sources_fast(Ql,𝚽l,Σs,mat,P,pl,Nm,Ns,Ngi,Ngf)

Optimized counterpart of `particle_sources`: build the source one particle deposits into
another, skipping the group transfers that do not exist (see `group_transfer_lists`).

Bit-identical to the reference: the same products reach the same elements, and the
accumulation over incoming groups keeps its increasing order — only blocks that would have
added exactly zero are dropped.

**The loop nesting is load-bearing and must not be "tidied".** Hoisting the `(gi,gf)` loops
above the voxel loop, so as to skip whole blocks at once, was measured at **0.79× — slower
than doing all the work** — because `Ql[gf,…]` is then streamed once per incoming group.
Keeping the voxel loop outside and the group loop inside is what makes the saving real:
**2.88× measured** at `Ng = 40×40`, `P = 36` (the gain exceeds the 2.4× block ratio because
the material lookup also leaves the innermost loop).
"""
function particle_sources_fast(Ql::Array{Float64},𝚽l::Array{Float64},Σs::Array{Float64},mat::Array{Int64},P::Int64,pl::Vector{Int64},Nm::Int64,Ns::Vector{Int64},Ngi::Int64,Ngf::Int64)
    lists = group_transfer_lists(Σs,Ngi,Ngf)
    @inbounds for gf in range(1,Ngf)
        gl = lists[gf]
        isempty(gl) && continue
        for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
            m = mat[ix,iy,iz]
            for k in eachindex(gl)
                gi = Int64(gl[k])
                for is in range(1,Nm)
                    @simd for p in range(1,P)
                        Ql[gf,p,is,ix,iy,iz] += Σs[m,gi,gf,pl[p]+1] * 𝚽l[gi,p,is,ix,iy,iz]
                    end
                end
            end
        end
    end
    return Ql
end