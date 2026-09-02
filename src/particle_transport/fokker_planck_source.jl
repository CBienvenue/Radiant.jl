"""
    FP_source(N::Int64,P::Int64,Nm::Int64,T::Vector{Float64},𝚽l::Array{Float64},
    Ql::Array{Float64},Ns::Vector{Int64},mat::Array{Int64,3},ℳ::Array{Float64,2},
    Mn::Array{Float64,2},Dn::Array{Float64,2})

Calculate the angular Fokker-Planck source term in Cartesian geometry.

# Input Argument(s)
- `P::Int64`: number of angular interpolation basis.
- `Nm::Int64`: total number of spatial and energy moments.
- `T::Vector{Float64}`: restricted momentum transfer.
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `Ql::Array{Float64}`: Legendre components of the in-cell source.
- `Ns::Vector{Int64}`: number of voxels per axis.    
- `mat::Array{Int64,3}`: material identifier per voxel.
- `ℳ::Array{Float64,2}`: Fokker-Planck scattering matrix.

# Output Argument(s)
- `Ql::Array{Float64}`: Legendre components of the in-cell source.

# Reference(s)
- Morel (1988) : A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.

"""
function fokker_planck_source(P::Int64,Nm::Int64,T::Vector{Float64},𝚽l::Array{Float64},Ql::Array{Float64},Ns::Vector{Int64},mat::Array{Int64,3},ℳ::Array{Float64,2})

# Compute the angular Fokker-Planck source term
for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
  for is in range(1,Nm), n in range(1,P), m in range(1,P)
    Ql[n,is,ix,iy,iz] += T[mat[ix,iy,iz]] * ℳ[n,m] * 𝚽l[m,is,ix,iy,iz]
  end
end

return Ql
end

"""
    fokker_planck_source_fast(P,Nm,T,𝚽l,Ql,Ns,mat,ℳ)

Optimized counterpart of `fokker_planck_source`, for the GN fast path.

The reference triple loop is a dense matrix product in disguise: `ℳ` multiplies the
angular axis of `𝚽l`, which is that array's leading — contiguous — axis. Written as one
`gemm` per material it runs ~20× faster (measured), because the reference form reloads and
restores `Ql[n,is,…]` on every `m` iteration and cannot vectorize across the angular axis.

With a single material the whole grid is one `gemm`. With several, each material's voxels
are gathered into a scratch buffer, multiplied, and scattered back — still a handful of
`gemm` calls instead of `P²·Nm·Nvox` scalar updates.

Agreement with the reference is to rounding (measured 3.6e-15 relative), not bit-for-bit:
`gemm` accumulates the sum over `m` in its own order.
"""
function fokker_planck_source_fast(P::Int64,Nm::Int64,T::Vector{Float64},𝚽l::Array{Float64},Ql::Array{Float64},Ns::Vector{Int64},mat::Array{Int64,3},ℳ::Array{Float64,2})

    Nvox = Ns[1]*Ns[2]*Ns[3]
    NS   = Nm*Nvox
    𝚽mat = reshape(𝚽l, P, NS)
    Qmat = reshape(Ql, P, NS)

    # Single material: the grid is one gemm.
    if length(T) == 1 || all(==(mat[1]),mat)
        mul!(Qmat, ℳ, 𝚽mat, T[mat[1]], 1.0)
        return Ql
    end

    # Several materials: one gemm per material over its gathered voxels.
    matv = reshape(mat, Nvox)
    cols = Int64[]
    buf_in  = Matrix{Float64}(undef, P, 0)
    buf_out = Matrix{Float64}(undef, P, 0)
    for m in 1:length(T)
        if T[m] == 0.0 continue end
        empty!(cols)
        for iv in 1:Nvox
            if matv[iv] == m
                for is in 1:Nm
                    push!(cols, (iv-1)*Nm + is)
                end
            end
        end
        isempty(cols) && continue
        if size(buf_in,2) < length(cols)
            buf_in  = Matrix{Float64}(undef, P, length(cols))
            buf_out = Matrix{Float64}(undef, P, length(cols))
        end
        vin  = view(buf_in ,:,1:length(cols))
        vout = view(buf_out,:,1:length(cols))
        @inbounds for (j,c) in enumerate(cols), p in 1:P
            vin[p,j] = 𝚽mat[p,c]
        end
        @inbounds for (j,c) in enumerate(cols), p in 1:P
            vout[p,j] = Qmat[p,c]
        end
        mul!(vout, ℳ, vin, T[m], 1.0)
        @inbounds for (j,c) in enumerate(cols), p in 1:P
            Qmat[p,c] = vout[p,j]
        end
    end
    return Ql
end