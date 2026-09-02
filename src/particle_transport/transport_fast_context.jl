"""
    transport_fast_context.jl

Precomputed context shared by the optimized solver chains (`gn_one_speed_fast` and below,
and their SN counterparts). It exists because everything the reference per-cell kernels
(`gn_3D_BTE!`, `gn_3D_BFP!`, `flux_3D_BTE`, …) recompute for every voxel is in fact
voxel-independent:

- the index maps (`index_Exyzp`, `index_Eyz`, …) are closures re-evaluated for every (i,j)
  pair of every cell — here they become flat `Int32` tables built once;
- the assembly coefficients (`C[i]/Δ`, `(-1)^n`, `s^(n-1)`, the `count(>(1),…)` activity
  tests) depend only on the scheme, the octant and the mesh widths;
- **the cell matrix `𝒮` itself** depends only on (Σt, S⁻, S⁺, S, Δx, Δy, Δz, C, ω, 𝒩,
  octant signs) — never on the voxel beyond its widths. Since `scheme_weights_gn` only
  admits DD and DG for GN, the closure weights `ω` are constant, so `𝒮` can be assembled
  and factorized once per (material, mesh-width triple, angular patch, energy group)
  instead of once per voxel: a few hundred factorizations instead of `Nvox × Npatch` per
  source-iteration pass.

Keeping LAPACK out of the per-cell loop matters on its own: `lu!`/`ldiv!` go through
OpenBLAS' globally-locked workspace allocator, so calling them per cell is far more expensive
than the arithmetic they perform.

The factorization itself still uses `lu!`, once per configuration and outside the hot
loop, so the factors are bit-identical to the reference chain's. Only the per-cell
triangular solves are hand-rolled, in the column-oriented order of reference BLAS `dtrsv`.

**Dimension- and solver-agnostic by construction.** The 3D moment tables reduce exactly to
the 2D ones at `Nmz = NmE = 1`, and to the 1D ones at `Nmy = Nmz = NmE = 1` — the reference
index maps `index_xy`/`index_xyp` of `gn_2D_*` and their 1D analogues are the 3D
`_gnf_index_*` evaluated at those degenerate sizes. Only the per-cell kernels and the sweeps
are dimension-specific, because a genuine 2D problem has no z streaming term at all (it is
not 3D with `Nz = 1`).

The SN correspondence is just as direct: an SN *direction* plays the role of a GN *patch*,
and an SN cell system is a GN one with `Nq = 1` — SN solves one direction at a time, so it
has no angular block. One caveat there: SN's adaptive scheme (`isAdapt`, AWD) recomputes the
closure weights *per cell* from the local flux, which makes the cell matrix voxel-dependent
and invalidates the factorization cache. The SN fast path must fall back to the reference
chain in that case. GN is unaffected — `scheme_weights_gn` only admits DD and DG.

The reference chains are untouched and remain the default; see `set_fast_path`.
"""

# ---------------------------------------------------------------------------------------
# In-cell moment tables
# ---------------------------------------------------------------------------------------

"""
    GNFastMoments

Flattened description of the in-cell moment space, shared by every angular patch and every
energy group. Replaces the `index_*` closures of the reference kernels.

`Nact` active moment tuples `(aE,ax,ay,az)` are listed in exactly the iteration order of
the reference kernels, so any accumulation performed in that order matches bit-for-bit.
For the fully-coupled case every tuple is active; otherwise a tuple is active when at most
one of its indices exceeds 1.

The system row of moment `k` and angular basis `ip` is `col[k] + (ip-1)*Nmom` — the
reference `index_Exyzp` is exactly `index_Exyz + (ip-1)*Nmom` in both the coupled and
uncoupled cases, so the unknown vector is the matrix `reshape(𝚽, Nmom, Nq)`.

The BTE tables are the BFP ones with `NmE = 1`: the reference `index_xy`/`index_yz`/
`index_xz`/`index_xyz` coincide with `index_Exy`/`index_Eyz`/`index_Exz`/`index_Exyz`
evaluated at `iE = 1, NmE = 1`.
"""
struct GNFastMoments{NMOM}
    Nmx  ::Int64
    Nmy  ::Int64
    Nmz  ::Int64
    NmE  ::Int64
    Nmom ::Int64                 # in-cell moments = columns of 𝚽n / rows of reshape(𝚽,·,Nq)
    Nact ::Int64                 # active moment tuples
    aE   ::Vector{Int32}         # (Nact) energy moment index of each active tuple
    ax   ::Vector{Int32}         # (Nact) x moment index
    ay   ::Vector{Int32}         # (Nact) y moment index
    az   ::Vector{Int32}         # (Nact) z moment index
    col  ::Vector{Int32}         # (Nact) index_Exyz  — column of 𝚽n / Qn
    cyz  ::Vector{Int32}         # (Nact) index_Eyz   — column of 𝚽x12
    cxz  ::Vector{Int32}         # (Nact) index_Exz   — column of 𝚽y12
    cxy  ::Vector{Int32}         # (Nact) index_Exy   — column of 𝚽z12
    cxyz ::Vector{Int32}         # (Nact) index_xyz   — column of 𝚽E12
end

# Index maps transcribed from gn_3D_bfp.jl / gn_3D_bte.jl, unchanged.
@inline function _gnf_index_Exy(iE,ix,iy,Nmx,NmE,isFC)
    if isFC
        return Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        return i
    end
end

@inline function _gnf_index_Exz(iE,ix,iz,Nmx,NmE,isFC)
    if isFC
        return Nmx*NmE*(iz-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iz > 1 i += NmE-1 + Nmx-1 end
        return i
    end
end

@inline function _gnf_index_Eyz(iE,iy,iz,Nmy,NmE,isFC)
    if isFC
        return Nmy*NmE*(iz-1) + NmE*(iy-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (iE-1)
        if iy > 1 i += NmE-1 end
        if iz > 1 i += NmE-1 + Nmy-1 end
        return i
    end
end

@inline function _gnf_index_xyz(ix,iy,iz,Nmx,Nmy,isFC)
    if isFC
        return Nmy*Nmx*(iz-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1)
        if iy > 1 i += Nmx-1 end
        if iz > 1 i += Nmx-1 + Nmy-1 end
        return i
    end
end

@inline function _gnf_index_Exyz(iE,ix,iy,iz,Nmx,Nmy,NmE,isFC)
    if isFC
        return Nmx*Nmy*NmE*(iz-1) + Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        if iz > 1 i += NmE-1 + Nmx-1 + Nmy-1 end
        return i
    end
end

"""
    gn_fast_moments(Nmx,Nmy,Nmz,NmE,isFC)

Build the in-cell moment tables. Pass `NmE = 1` for the BTE kernels.

The tuples are enumerated in the reference kernels' order — `for ix, iy, iz, iE` with the
inactive ones skipped — so downstream accumulations reproduce the reference summation order.
"""
function gn_fast_moments(Nmx::Int64,Nmy::Int64,Nmz::Int64,NmE::Int64,isFC::Bool)

    Nmom = isFC ? Nmx*Nmy*Nmz*NmE : (Nmx+Nmy+Nmz+NmE-3)

    aE = Int32[]; ax = Int32[]; ay = Int32[]; az = Int32[]
    col = Int32[]; cyz = Int32[]; cxz = Int32[]; cxy = Int32[]; cxyz = Int32[]

    for ix in 1:Nmx, iy in 1:Nmy, iz in 1:Nmz, iE in 1:NmE
        if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
        push!(aE,iE); push!(ax,ix); push!(ay,iy); push!(az,iz)
        push!(col , _gnf_index_Exyz(iE,ix,iy,iz,Nmx,Nmy,NmE,isFC))
        push!(cyz , _gnf_index_Eyz(iE,iy,iz,Nmy,NmE,isFC))
        push!(cxz , _gnf_index_Exz(iE,ix,iz,Nmx,NmE,isFC))
        push!(cxy , _gnf_index_Exy(iE,ix,iy,Nmx,NmE,isFC))
        push!(cxyz, _gnf_index_xyz(ix,iy,iz,Nmx,Nmy,isFC))
    end

    # Nmom is carried in the type so the per-voxel kernels see the moment count as a
    # compile-time constant and unroll their loops over it (the index tables below are the
    # same length: `col` is a bijection onto 1:Nmom).
    return GNFastMoments{Nmom}(Nmx,Nmy,Nmz,NmE,Nmom,length(aE),aE,ax,ay,az,col,cyz,cxz,cxy,cxyz)
end

# ---------------------------------------------------------------------------------------
# Closure relations
# ---------------------------------------------------------------------------------------

"""
    GNFastClosure

One outgoing-face (or outgoing-energy) closure relation, flattened into CSR-like form.

The reference kernels write, for each face column `face[k]`,

    Φface[ip,face[k]] = ω0 * Φface[ip,face[k]]
                      + Σ over the entries of slice k of  coef * Φmat[row,ip]

with the entries enumerated along the closed-over axis in increasing order — the order
kept here, so the accumulation is bit-identical.
"""
struct GNFastClosure
    ω0    ::Vector{Float64}      # (nface) scaling of the incoming face value
    face  ::Vector{Int32}        # (nface) column in the face array
    off   ::Vector{Int32}        # (nface+1) offsets into row/coef
    row   ::Vector{Int32}        # moment index in reshape(𝚽,Nmom,Nq)
    coef  ::Vector{Float64}
end

# Build a closure by sweeping the three "kept" indices and accumulating over the fourth.
# `keep` selects which axis is closed over: :E, :x, :y or :z.
function _gn_fast_closure(keep::Symbol,mom::GNFastMoments,isFC::Bool,ω::Vector{Float64},
                          C::Vector{Float64},s::Int64)

    Nmx,Nmy,Nmz,NmE = mom.Nmx,mom.Nmy,mom.Nmz,mom.NmE
    face = Int32[]; off = Int32[1]; row = Int32[]; coef = Float64[]

    # Outer tuples in the reference kernels' loop order (`for a, b, c` nests with `a`
    # slowest, which is what the nested-`for` comprehension below reproduces), and the
    # length of the axis being accumulated over.
    outer, naxis = if keep === :E
        ([(ix,iy,iz) for ix in 1:Nmx for iy in 1:Nmy for iz in 1:Nmz], NmE)
    elseif keep === :x
        ([(iE,iy,iz) for iE in 1:NmE for iy in 1:Nmy for iz in 1:Nmz], Nmx)
    elseif keep === :y
        ([(iE,ix,iz) for iE in 1:NmE for ix in 1:Nmx for iz in 1:Nmz], Nmy)
    else
        ([(iE,ix,iy) for iE in 1:NmE for ix in 1:Nmx for iy in 1:Nmy], Nmz)
    end

    for t in outer
        # The reference skips an outer tuple when two of its indices exceed 1.
        if (~isFC) && (count(>(1),t) ≥ 2) continue end

        # Face column and the full 4-tuple as a function of the accumulation index.
        fcol, tuple_of = if keep === :E
            (ix,iy,iz) = t
            (_gnf_index_xyz(ix,iy,iz,Nmx,Nmy,isFC), n -> (n,ix,iy,iz))
        elseif keep === :x
            (iE,iy,iz) = t
            (_gnf_index_Eyz(iE,iy,iz,Nmy,NmE,isFC), n -> (iE,n,iy,iz))
        elseif keep === :y
            (iE,ix,iz) = t
            (_gnf_index_Exz(iE,ix,iz,Nmx,NmE,isFC), n -> (iE,ix,n,iz))
        else
            (iE,ix,iy) = t
            (_gnf_index_Exy(iE,ix,iy,Nmx,NmE,isFC), n -> (iE,ix,iy,n))
        end

        push!(face, fcol)
        for n in 1:naxis
            (iE,ix,iy,iz) = tuple_of(n)
            if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
            push!(row, _gnf_index_Exyz(iE,ix,iy,iz,Nmx,Nmy,NmE,isFC))
            # Coefficients transcribed from the reference closure relations.
            c = if keep === :E
                C[n] * (-1)^(n-1) * ω[n+1]
            else
                C[n] * s^(n-1) * ω[n+1]
            end
            push!(coef, c)
        end
        push!(off, Int32(length(row)+1))
    end

    return GNFastClosure(fill(ω[1],length(face)),face,off,row,coef)
end

# ---------------------------------------------------------------------------------------
# Per-cell system cache
# ---------------------------------------------------------------------------------------

"""
    GNFastCells

Maps a voxel to the index of its cell-system configuration. A configuration is a
(material, Δx, Δy, Δz) combination; the per-cell matrix depends on nothing else.

`Nconf` stays small on the meshes Radiant builds (uniform within each region), which is
what makes the factorization cache worthwhile. `gn_fast_context` refuses to build a
context when it would not.
"""
struct GNFastCells
    Nconf ::Int64
    Nmat  ::Int64
    Nkx   ::Int64
    Nky   ::Int64
    Nkz   ::Int64
    kx    ::Vector{Int32}        # (Nx) distinct-width index of each x slab
    ky    ::Vector{Int32}        # (Ny)
    kz    ::Vector{Int32}        # (Nz)
    Δx    ::Vector{Float64}      # (Nkx) the distinct widths themselves
    Δy    ::Vector{Float64}
    Δz    ::Vector{Float64}
end

# conf = mat + Nmat*((kx-1) + Nkx*((ky-1) + Nky*(kz-1)))
@inline function gn_fast_conf(cells::GNFastCells,m::Int64,ix::Int64,iy::Int64,iz::Int64)
    @inbounds return m + cells.Nmat*((cells.kx[ix]-1) +
                        cells.Nkx*((cells.ky[iy]-1) + cells.Nky*(cells.kz[iz]-1)))
end

function _gn_fast_distinct(Δ::Vector{Float64})
    vals = Float64[]
    idx  = Vector{Int32}(undef,length(Δ))
    for (i,d) in enumerate(Δ)
        k = findfirst(==(d),vals)
        if k === nothing
            push!(vals,d); k = length(vals)
        end
        idx[i] = Int32(k)
    end
    return vals, idx
end

# ---------------------------------------------------------------------------------------
# Per-patch context
# ---------------------------------------------------------------------------------------

"""
    GNFastPatch

Everything one angular patch (octant `u`, subdivision `v,w`) needs, for one energy group:
the LU factors of every cell-system configuration, and the flattened assembly and closure
coefficients.
"""
struct GNFastPatch{NQ}
    Nq     ::Int64
    Nm     ::Int64               # system size = Nmom * Nq
    sx     ::Int64
    sy     ::Int64
    sz     ::Int64
    𝒩x     ::Matrix{Float64}
    𝒩y     ::Matrix{Float64}
    𝒩z     ::Matrix{Float64}
    LU     ::Array{Float64,3}    # (Nm,Nm,Nconf) factors from lu!
    ipiv   ::Matrix{Int64}       # (Nm,Nconf)
    qcx    ::Matrix{Float64}     # (Nact,Nkx) incoming-x coefficient per active moment
    qcy    ::Matrix{Float64}     # (Nact,Nky)
    qcz    ::Matrix{Float64}     # (Nact,Nkz)
    qcE    ::Matrix{Float64}     # (Nact,Nmat) incoming-energy coefficient (BFP; empty for BTE)
    clE    ::GNFastClosure
    clx    ::GNFastClosure
    cly    ::GNFastClosure
    clz    ::GNFastClosure
end

"""
    GNFastContext

The whole optimized-chain context for one energy group: moment tables, voxel→configuration
map, and one `GNFastPatch` per active angular patch (indexed by the linear patch id used by
`gn_inner_pass_fast!`).
"""
struct GNFastContext{NMOM,NQ}
    mom    ::GNFastMoments{NMOM}
    cells  ::GNFastCells
    patch  ::Vector{GNFastPatch{NQ}}
    isCSD  ::Bool
end

# ---------------------------------------------------------------------------------------
# Per-sweep workspace
# ---------------------------------------------------------------------------------------

"""
    GNFastWorkspace

Scratch for one patch sweep. The reference chain keeps a single shared set of these
(`gn_one_speed.jl`); the optimized chain allocates one per patch so a sweep never has to
clear state left by the previous one.
"""
struct GNFastWorkspace
    Q    ::Vector{Float64}
    𝚽    ::Vector{Float64}
    bufx ::Vector{Float64}       # moving x-boundary, shape (Nq, Nmf[1]) — innermost sweep axis
    bufy ::Vector{Float64}       # moving y-boundary, shape (Nq, Nmf[2], Nx)
    bufz ::Vector{Float64}       # moving z-boundary, shape (Nq, Nmf[3], Nx, Ny)
    _pad ::Vector{Float64}
end

"""
    GNFastWorkspace(Ndims,Nm,Nq,Nmf,Ns)

Allocate one sweep workspace. The moving buffers are sized by their role in the sweep
nesting, which puts the fastest array index innermost: the innermost axis carries a single
face, the next one a line, the outermost a plane. Axes the geometry does not have get no
buffer at all.
"""
function GNFastWorkspace(Ndims::Int64,Nm::Int64,Nq::Int64,Nmf::Vector{Int64},Ns::Vector{Int64})
    nx = Ndims ≥ 1 ? Nq*Nmf[1] : 0                       # innermost sweep axis: one face
    ny = Ndims ≥ 2 ? Nq*Nmf[2]*Ns[1] : 0                 # middle axis: a line along x
    nz = Ndims ≥ 3 ? Nq*Nmf[3]*Ns[1]*Ns[2] : 0           # outer axis: an x-y plane
    return GNFastWorkspace(zeros(Nm),zeros(Nm),zeros(nx),zeros(ny),zeros(nz),zeros(8))
end

# ---------------------------------------------------------------------------------------
# Context construction
# ---------------------------------------------------------------------------------------

"""
    gn_fast_applicable(Ndims,Δs,mat,Nmat,𝒪,Nq,isFC;max_cache_bytes)

Whether the optimized chain can serve this problem, and why not when it cannot.

The factorization cache is only worthwhile when the number of distinct
(material, Δx, Δy, Δz) configurations stays small — true on the region-uniform meshes
`Geometry` builds, and on the uniform grids used for dose calculations. Returns
`(true,"")` or `(false,reason)`.
"""
function gn_fast_applicable(Ndims::Int64,Δs::Vector{Vector{Float64}},Nmat::Int64,
                            𝒪::Vector{Int64},Nq::Int64,isFC::Bool,Npatch::Int64;
                            max_cache_bytes::Int64=1<<29)

    if Ndims ∉ (1,2,3)
        return false, "the optimized chain covers 1D, 2D and 3D"
    end

    Nkx = length(unique(Δs[1]))
    Nky = Ndims ≥ 2 ? length(unique(Δs[2])) : 1
    Nkz = Ndims ≥ 3 ? length(unique(Δs[3])) : 1
    Nconf = Nmat*Nkx*Nky*Nkz
    Nmom  = isFC ? prod(𝒪) : (1+sum(𝒪.-1))
    Nm    = Nmom*Nq
    bytes = Npatch*Nconf*Nm*Nm*8

    if bytes > max_cache_bytes
        return false, "the cell-system cache would need " *
                      "$(round(bytes/2^30,digits=2)) GiB " *
                      "($Nconf material/mesh-width configurations × $Npatch patches × " *
                      "$(Nm)² ) — use a mesh that is uniform within each region"
    end
    return true, ""
end

"""
    gn_fast_context(...)

Build the optimized-chain context for one energy group.

`Σt`, `S⁻`, `S⁺`, `S` are the per-material group constants; `𝒩` holds the patch streaming
weight matrices, `patch_uvw` the (u,v,w) triple of each linear patch id. The matrices are
assembled by `gn_3D_BFP_matrix!` / `gn_3D_BTE_matrix!`, which carry the reference assembly
loops verbatim, then factorized with `lu!` — the same factorization the reference chain
performs per cell.
"""
function gn_fast_context(Ndims::Int64,Δs::Vector{Vector{Float64}},
                         Nmat::Int64,Σt::Vector{Float64},S⁻::Vector{Float64},
                         S⁺::Vector{Float64},S::Array{Float64},
                         𝒪::Vector{Int64},isFC::Bool,isCSD::Bool,
                         C::Vector{Float64},ω::Vector{Vector{Float64}},𝒲::Array{Float64},
                         Nq::Int64,𝒩::Array{Float64},
                         patch_uvw::Vector{NTuple{3,Int64}})

    # A genuine 2D (resp. 1D) problem has no z (resp. no y,z) streaming term, so the unused
    # axes are collapsed to a single moment. `get_schemes` already sets 𝒪 = 1 on them; this
    # makes the intent explicit and keeps the degenerate closures out of the context.
    Nmx = 𝒪[1]
    Nmy = Ndims ≥ 2 ? 𝒪[2] : 1
    Nmz = Ndims ≥ 3 ? 𝒪[3] : 1
    NmE = isCSD ? 𝒪[4] : 1
    mom = gn_fast_moments(Nmx,Nmy,Nmz,NmE,isFC)
    Nm  = mom.Nmom*Nq

    # Distinct mesh widths per axis
    # `Δs` only carries the axes the geometry has; the others are undefined references. An
    # inactive axis has a single, irrelevant width — no kernel streams along it — so a
    # placeholder keeps the configuration map correct and minimal.
    Δxv,kx = _gn_fast_distinct(Δs[1])
    Δyv,ky = Ndims ≥ 2 ? _gn_fast_distinct(Δs[2]) : ([1.0],Int32[1])
    Δzv,kz = Ndims ≥ 3 ? _gn_fast_distinct(Δs[3]) : ([1.0],Int32[1])
    Nkx,Nky,Nkz = length(Δxv),length(Δyv),length(Δzv)
    Nconf = Nmat*Nkx*Nky*Nkz
    cells = GNFastCells(Nconf,Nmat,Nkx,Nky,Nkz,kx,ky,kz,Δxv,Δyv,Δzv)

    sxs = [1,1,1,1,-1,-1,-1,-1]
    sys = [1,1,-1,-1,1,1,-1,-1]
    szs = [1,-1,1,-1,1,-1,1,-1]

    # g(n,s) of the reference kernels
    gsign(n,s) = s > 0 ? 1 : -(-1)^(n-1)

    patches = Vector{GNFastPatch{Nq}}(undef,length(patch_uvw))

    for (ip_lin,(u,v,w)) in enumerate(patch_uvw)
        sx,sy,sz = sxs[u],sys[u],szs[u]
        # `𝒩` carries one streaming-weight matrix per axis the geometry has, so the z slot
        # simply does not exist in 2D.
        𝒩x = Matrix{Float64}(𝒩[:,:,1,u,v,w])
        𝒩y = Ndims ≥ 2 ? Matrix{Float64}(𝒩[:,:,2,u,v,w]) : zeros(Nq,Nq)
        𝒩z = Ndims ≥ 3 ? Matrix{Float64}(𝒩[:,:,3,u,v,w]) : zeros(Nq,Nq)

        LU   = Array{Float64,3}(undef,Nm,Nm,Nconf)
        ipiv = Matrix{Int64}(undef,Nm,Nconf)
        𝒮    = zeros(Nm,Nm)

        for m in 1:Nmat, ikx in 1:Nkx, iky in 1:Nky, ikz in 1:Nkz
            conf = m + Nmat*((ikx-1) + Nkx*((iky-1) + Nky*(ikz-1)))
            if Ndims == 3
                if isCSD
                    gn_3D_BFP_matrix!(𝒮,sx,sy,sz,Σt[m],S⁺[m],view(S,m,:),
                                      Δxv[ikx],Δyv[iky],Δzv[ikz],
                                      Nmx,Nmy,Nmz,NmE,Nq,C,ω[1],ω[2],ω[3],ω[4],
                                      𝒩x,𝒩y,𝒩z,𝒲,isFC)
                else
                    gn_3D_BTE_matrix!(𝒮,sx,sy,sz,Σt[m],Δxv[ikx],Δyv[iky],Δzv[ikz],
                                      Nmx,Nmy,Nmz,Nq,C,ω[1],ω[2],ω[3],𝒩x,𝒩y,𝒩z,isFC)
                end
            elseif Ndims == 2
                if isCSD
                    gn_2D_BFP_matrix!(𝒮,sx,sy,Σt[m],S⁺[m],view(S,m,:),
                                      Δxv[ikx],Δyv[iky],
                                      Nmx,Nmy,NmE,Nq,C,ω[1],ω[2],ω[4],𝒩x,𝒩y,𝒲,isFC)
                else
                    gn_2D_BTE_matrix!(𝒮,sx,sy,Σt[m],Δxv[ikx],Δyv[iky],
                                      Nmx,Nmy,Nq,C,ω[1],ω[2],𝒩x,𝒩y,isFC)
                end
            elseif Ndims == 1
                if isCSD
                    gn_1D_BFP_matrix!(𝒮,sx,Σt[m],S⁺[m],view(S,m,:),Δxv[ikx],
                                      Nmx,NmE,Nq,C,ω[1],ω[4],𝒩x,𝒲,isFC)
                else
                    gn_1D_BTE_matrix!(𝒮,sx,Σt[m],Δxv[ikx],Nmx,Nq,C,ω[1],𝒩x)
                end
            else
                error("gn_fast_context: unsupported dimension $(Ndims).")
            end
            F = lu!(𝒮)
            LU[:,:,conf]  = F.factors
            ipiv[:,conf] .= F.ipiv
        end

        # Right-hand-side coefficients, transcribed from the reference source loops.
        qcx = zeros(mom.Nact,Nkx); qcy = zeros(mom.Nact,Nky); qcz = zeros(mom.Nact,Nkz)
        qcE = isCSD ? zeros(mom.Nact,Nmat) : zeros(0,0)
        for k in 1:mom.Nact
            ix,iy,iz,iE = mom.ax[k],mom.ay[k],mom.az[k],mom.aE[k]
            for ikx in 1:Nkx
                qcx[k,ikx] = -C[ix]/Δxv[ikx] * (gsign(ix,sx)*ω[1][1] + gsign(ix,-sx))
            end
            for iky in 1:Nky
                qcy[k,iky] = -C[iy]/Δyv[iky] * (gsign(iy,sy)*ω[2][1] + gsign(iy,-sy))
            end
            for ikz in 1:Nkz
                qcz[k,ikz] = -C[iz]/Δzv[ikz] * (gsign(iz,sz)*ω[3][1] + gsign(iz,-sz))
            end
            if isCSD
                for m in 1:Nmat
                    qcE[k,m] = C[iE] * ((-1)^iE*S⁺[m]*ω[4][1] + S⁻[m])
                end
            end
        end

        # Closure relations. The energy closure is only meaningful for the BFP kernel, and
        # the y/z ones only for the dimensions that actually stream along those axes; the
        # others are built empty so the kernels can skip them without a branch.
        ωE = isCSD ? ω[4] : [0.0,1.0]
        empty_cl = GNFastClosure(Float64[],Int32[],Int32[1],Int32[],Float64[])
        clE = _gn_fast_closure(:E,mom,isFC,ωE,C,1)
        clx = _gn_fast_closure(:x,mom,isFC,ω[1],C,sx)
        cly = Ndims ≥ 2 ? _gn_fast_closure(:y,mom,isFC,ω[2],C,sy) : empty_cl
        clz = Ndims ≥ 3 ? _gn_fast_closure(:z,mom,isFC,ω[3],C,sz) : empty_cl

        patches[ip_lin] = GNFastPatch{Nq}(Nq,Nm,sx,sy,sz,𝒩x,𝒩y,𝒩z,LU,ipiv,
                                          qcx,qcy,qcz,qcE,clE,clx,cly,clz)
    end

    return GNFastContext{mom.Nmom,Nq}(mom,cells,patches,isCSD)
end

# ---------------------------------------------------------------------------------------
# Per-cell triangular solve
# ---------------------------------------------------------------------------------------

"""
    gn_fast_chunks(n,nchunk)

Split `1:n` into at most `nchunk` contiguous ranges of near-equal length.

Used to cut the angular transforms along the spatial axis, so that every patch's transform
for a given range touches the same slice of the full-range array: it then stays in cache
across the whole patch loop instead of being streamed from memory once per patch. The chain
being serial, a single range is the usual case; the split is kept for cache blocking.
"""
function gn_fast_chunks(n::Int64,nchunk::Int64)
    nchunk = max(1,min(nchunk,n))
    q, r = divrem(n,nchunk)
    chunks = Vector{UnitRange{Int64}}(undef,nchunk)
    lo = 1
    for i in 1:nchunk
        len = q + (i ≤ r ? 1 : 0)
        chunks[i] = lo:(lo+len-1)
        lo += len
    end
    return chunks
end

"""
    GN_FAST_BLOCK

Number of spatial columns processed at a time by the batched angular transforms
(`gn_fast_scatter!` / `gn_fast_gather!`). It sizes their staging buffer at
`Nq·Npatch × GN_FAST_BLOCK`, which has to stay small enough to live in cache — the whole
point of blocking is that `Nq·Npatch × NS` would be 8.6 GB on a 100³ mesh.

Measured flat between 512 and 5000 columns (2.57× to 2.65× against the unbatched form), so
the value is chosen for the buffer size rather than for speed: 1024 columns give ~1.7 MB.
"""
const GN_FAST_BLOCK = 1024

"""
    gn_fast_scatter!(dst,tmp,rng,Nq,Npatch,NS)

Scatter one block of the batched forward transform into the per-patch layout.

The batched `gemm` produces `tmp[(k-1)·Nq+q, c]` — patch-major — while the sweep needs
`[q, cell, patch]`, where a patch's slice is contiguous. This is the copy that reconciles
the two; it is the price of keeping the sweep's layout, and the 2.59× measured for the
batched transform already includes it.
"""
@inline function gn_fast_scatter!(dst::Vector{Float64},tmp::Matrix{Float64},
                                  rng::UnitRange{Int64},Nq::Int64,Npatch::Int64,NS::Int64)
    lo = first(rng); nb = length(rng)
    @inbounds for k in 1:Npatch
        bd = Nq*NS*(k-1) + Nq*(lo-1)
        bt = Nq*(k-1)
        for c in 1:nb
            d = bd + Nq*(c-1)
            for q in 1:Nq
                dst[d+q] = tmp[bt+q,c]
            end
        end
    end
    return nothing
end

"""
    gn_fast_gather!(tmp,src,rng,Nq,Npatch,NS)

Gather one block of the per-patch layout into the patch-major staging buffer, so the inverse
transform can also be a single `gemm`. The mirror of `gn_fast_scatter!`.
"""
@inline function gn_fast_gather!(tmp::Matrix{Float64},src::Vector{Float64},
                                 rng::UnitRange{Int64},Nq::Int64,Npatch::Int64,NS::Int64)
    lo = first(rng); nb = length(rng)
    @inbounds for k in 1:Npatch
        bs = Nq*NS*(k-1) + Nq*(lo-1)
        bt = Nq*(k-1)
        for c in 1:nb
            s = bs + Nq*(c-1)
            for q in 1:Nq
                tmp[bt+q,c] = src[s+q]
            end
        end
    end
    return nothing
end

"""
    gn_fast_closure!(Φ,off,cl,𝚽,Nmom,Nq)

Apply one flattened closure relation (see `GNFastClosure`) to an outgoing face — or to the
outgoing energy flux — of a voxel.

`Φ` is the flat storage of the face array and `off` the voxel's base offset in it; the face
array has the angular index leading, so element `(ip, c)` sits at `off + ip + Nq*(c-1)`.
Taking a flat vector and an offset rather than a 2-D slice keeps the sweep from building a
`SubArray` per voxel — measured at ~16 % of the whole run before this change.

The reference kernels run `for ip in 1:Np` outside all closure loops; here the angular index
is innermost. Each `(ip, face column)` accumulation is independent of every other, and the
entries of a slice are stored in the reference order, so the summation each output element
sees is unchanged.
"""
@inline function gn_fast_closure!(Φ::Vector{Float64},off::Int64,cl::GNFastClosure,
                                  𝚽::Vector{Float64},::Val{NMOM},::Val{NQ}) where {NMOM,NQ}
    Nmom = NMOM; Nq = NQ
    @inbounds for kf in 1:length(cl.face)
        b  = off + Nq*(cl.face[kf]-1)
        ω0 = cl.ω0[kf]
        for ip in 1:Nq
            Φ[b+ip] *= ω0
        end
        for e in cl.off[kf]:cl.off[kf+1]-1
            r = cl.row[e]; a = cl.coef[e]
            for ip in 1:Nq
                Φ[b+ip] += a * 𝚽[r + (ip-1)*Nmom]
            end
        end
    end
    return nothing
end

"""
    gn_fast_solve!(𝚽,LU,ipiv,Q,Nm,conf)

Solve `A x = b` from the cached factors, replacing the reference chain's
`ldiv!(𝚽,lu!(𝒮),Q)` (LAPACK `dgetrs`) without entering LAPACK — whose locked workspace
allocator costs more per call than the arithmetic on systems this small.

Reproduces `dgetrs('N',…)` step for step: row interchanges, unit-lower forward
substitution, then upper back substitution, both in the column-oriented order of reference
BLAS `dtrsv` for column-major storage.
"""
@inline function gn_fast_solve!(𝚽::Vector{Float64},LU::Array{Float64,3},
                                ipiv::Matrix{Int64},Q::Vector{Float64},
                                Nm::Int64,conf::Int64)
    @inbounds begin
        # Row interchanges (dlaswp)
        for k in 1:Nm
            p = ipiv[k,conf]
            if p != k
                t = Q[k]; Q[k] = Q[p]; Q[p] = t
            end
        end
        # Forward substitution, unit lower triangle (dtrsv 'L','N','U')
        for j in 1:Nm
            xj = Q[j]
            if xj != 0.0
                for i in j+1:Nm
                    Q[i] -= xj * LU[i,j,conf]
                end
            end
        end
        # Back substitution, upper triangle (dtrsv 'U','N','N')
        for j in Nm:-1:1
            xj = Q[j] / LU[j,j,conf]
            Q[j] = xj
            if xj != 0.0
                for i in 1:j-1
                    Q[i] -= xj * LU[i,j,conf]
                end
            end
        end
        for i in 1:Nm
            𝚽[i] = Q[i]
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# SN context
# ---------------------------------------------------------------------------------------

"""
    _sn_fast_closure(keep,mom,isFC,ω,C,s)

Flatten one SN closure relation, the discrete-ordinates counterpart of `_gn_fast_closure`.

Two things differ from GN. SN's closure weights carry three indices — `ωx[jx+1,jy,jz]` —
because the WD/AWD family couples the transverse moments, so the incoming-face scaling is
per face column rather than a single constant. And SN solves one direction at a time, so the
system index is the moment index itself: there is no angular block.

The same routine serves the energy axis of the BFP kernel (`keep = :E`), whose closure has
the same shape with `s = -1`; the weights then carry four indices instead of three.

The reference applies the scaling inside the `for jx, jy, jz` nest, guarded by `jx == 1`;
since `jx` is the outermost loop that is the first visit of each face column, and the
accumulation that follows runs over increasing `jx` — the order kept here.
"""
function _sn_fast_closure(keep::Symbol,mom::GNFastMoments,isFC::Bool,ω::Array{Float64},
                          C::Vector{Float64},s::Float64)

    Nmx,Nmy,Nmz,NmE = mom.Nmx,mom.Nmy,mom.Nmz,mom.NmE
    # BTE weights carry three indices, BFP weights four; give both the same shape so the
    # lookups below are written once. In BTE `NmE == 1`, so the added index is always 1.
    ωr = ndims(ω) == 3 ? reshape(ω,size(ω,1),size(ω,2),size(ω,3),1) : ω
    ω0 = Float64[]; face = Int32[]; off = Int32[1]; row = Int32[]; coef = Float64[]

    outer, naxis = if keep === :x
        ([(jy,jz,jE) for jy in 1:Nmy for jz in 1:Nmz for jE in 1:NmE], Nmx)
    elseif keep === :y
        ([(jx,jz,jE) for jx in 1:Nmx for jz in 1:Nmz for jE in 1:NmE], Nmy)
    elseif keep === :z
        ([(jx,jy,jE) for jx in 1:Nmx for jy in 1:Nmy for jE in 1:NmE], Nmz)
    else
        ([(jx,jy,jz) for jx in 1:Nmx for jy in 1:Nmy for jz in 1:Nmz], NmE)
    end

    for t in outer
        # A face column exists as soon as one moment tuple through it is active.
        if (~isFC) && (count(>(1),t) ≥ 2) continue end
        fcol, tuple_of, ω0v = if keep === :x
            (jy,jz,jE) = t
            (_gnf_index_Eyz(jE,jy,jz,Nmy,NmE,isFC), n -> (n,jy,jz,jE), ωr[1,jy,jz,jE])
        elseif keep === :y
            (jx,jz,jE) = t
            (_gnf_index_Exz(jE,jx,jz,Nmx,NmE,isFC), n -> (jx,n,jz,jE), ωr[1,jx,jz,jE])
        elseif keep === :z
            (jx,jy,jE) = t
            (_gnf_index_Exy(jE,jx,jy,Nmx,NmE,isFC), n -> (jx,jy,n,jE), ωr[1,jx,jy,jE])
        else
            (jx,jy,jz) = t
            (_gnf_index_xyz(jx,jy,jz,Nmx,Nmy,isFC), n -> (jx,jy,jz,n), ωr[1,jx,jy,jz])
        end

        push!(face,fcol); push!(ω0,ω0v)
        for n in 1:naxis
            (jx,jy,jz,jE) = tuple_of(n)
            if (~isFC) && (count(>(1),(jE,jx,jy,jz)) ≥ 2) continue end
            push!(row, _gnf_index_Exyz(jE,jx,jy,jz,Nmx,Nmy,NmE,isFC))
            c = if keep === :x
                C[jx] * s^(jx-1) * ωr[jx+1,jy,jz,jE]
            elseif keep === :y
                C[jy] * s^(jy-1) * ωr[jy+1,jx,jz,jE]
            elseif keep === :z
                C[jz] * s^(jz-1) * ωr[jz+1,jx,jy,jE]
            else
                C[jE] * s^(jE-1) * ωr[jE+1,jx,jy,jz]
            end
            push!(coef,c)
        end
        push!(off, Int32(length(row)+1))
    end

    return GNFastClosure(ω0,face,off,row,coef)
end

"""
    SNFastDir

Everything one discrete ordinate needs, for one energy group: the LU factors of every
cell-system configuration, and the flattened right-hand-side and closure coefficients.

The SN counterpart of `GNFastPatch`. The cell system is a GN one with `Nq = 1`, so its size
is the moment count alone.
"""
struct SNFastDir{NMOM}
    Nm   ::Int64
    sx   ::Int64                 # sweep ordering signs — `μ ≥ 0 ? 1 : -1`, as the reference
                                 # sweep tests it, which differs from `sign(μ)` at μ = 0
    sy   ::Int64
    sz   ::Int64
    LU   ::Array{Float64,3}      # (Nm,Nm,Nconf)
    ipiv ::Matrix{Int64}
    qcx  ::Matrix{Float64}       # (NMOM,Nkx) incoming-x coefficient per moment tuple
    qcy  ::Matrix{Float64}       # (NMOM,Nky)
    qcz  ::Matrix{Float64}       # (NMOM,Nkz)
    qcE  ::Matrix{Float64}       # (NMOM,Nmat) incoming-energy coefficient (BFP only)
    clx  ::GNFastClosure
    cly  ::GNFastClosure
    clz  ::GNFastClosure
    clE  ::GNFastClosure         # energy-axis closure (BFP only; empty for BTE)
end

"""
    SNFastContext

The optimized SN chain's context for one energy group: the moment tables, the
voxel→configuration map, and one `SNFastDir` per discrete ordinate.
"""
struct SNFastContext{NMOM}
    mom  ::GNFastMoments{NMOM}
    cells::GNFastCells
    dir  ::Vector{SNFastDir{NMOM}}
end

"""
    _sn_closure_2D_bte(keep,mom,isFC,ω,C,s)

Flatten one SN closure relation of the 2D BTE kernel.

The 2D BTE kernel does not have the same closure shape as the others: its weights carry the
*row's* transverse index as well as the column's — `ωx[jx+1,jy,iy]` — and the closure sums over
it. That is the full transverse coupling of the WD family; the 3D kernels and the 2D BFP one
keep only the column index. The `𝒪y` (resp. `𝒪x`) contributions are pushed as separate CSR
entries on the same row rather than summed first, so the accumulation order is the reference's.
"""
function _sn_closure_2D_bte(keep::Symbol,mom::GNFastMoments,isFC::Bool,ω::Array{Float64},
                            C::Vector{Float64},s::Float64)

    Nmx,Nmy = mom.Nmx,mom.Nmy
    ω0 = Float64[]; face = Int32[]; off = Int32[1]; row = Int32[]; coef = Float64[]
    outer, naxis, ncross = keep === :x ? (1:Nmy,Nmx,Nmy) : (1:Nmx,Nmy,Nmx)

    for t in outer
        fcol = keep === :x ? _gnf_index_Eyz(1,t,1,Nmy,1,isFC) : _gnf_index_Exz(1,t,1,Nmx,1,isFC)
        push!(face,fcol); push!(ω0,ω[1,t,t])
        for n in 1:naxis
            jx,jy = keep === :x ? (n,t) : (t,n)
            if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
            r = _gnf_index_Exyz(1,jx,jy,1,Nmx,Nmy,1,isFC)
            for i in 1:ncross
                push!(row,r)
                push!(coef, C[n] * s^(n-1) * ω[n+1,t,i])
            end
        end
        push!(off, Int32(length(row)+1))
    end

    return GNFastClosure(ω0,face,off,row,coef)
end

"""
    _sn_closure_1D(keep,mom,isFC,ω,C,s,isTangential)

Flatten one SN closure relation of the 1D kernels.

`flux_1D_BTE` has scalar weights and a scalar face flux. `flux_1D_BFP` has the 2D BTE shape on
both of its axes — `ωx[jx+1,jE,iE]` and `ωE[jE+1,jx,ix]`, each summed over the row index. A
tangential direction (`|μ| ≤ 1e-10`) skips the spatial closure entirely, as the reference does,
and gets an empty relation here.
"""
function _sn_closure_1D(keep::Symbol,mom::GNFastMoments,isFC::Bool,ω::Array{Float64},
                        C::Vector{Float64},s::Float64,isTangential::Bool)

    Nmx,NmE = mom.Nmx,mom.NmE
    ω0 = Float64[]; face = Int32[]; off = Int32[1]; row = Int32[]; coef = Float64[]
    if keep === :x && isTangential
        return GNFastClosure(ω0,face,off,row,coef)
    end
    scalar = ndims(ω) == 1                       # BTE: a single weight vector
    outer, naxis, ncross = keep === :x ? (1:NmE,Nmx,NmE) : (1:Nmx,NmE,Nmx)

    for t in outer
        fcol = keep === :x ? _gnf_index_Eyz(t,1,1,1,NmE,isFC) : _gnf_index_xyz(t,1,1,Nmx,1,isFC)
        push!(face,fcol); push!(ω0, scalar ? ω[1] : ω[1,t,t])
        for n in 1:naxis
            jx,jE = keep === :x ? (n,t) : (t,n)
            if (~isFC) && (count(>(1),(jE,jx)) ≥ 2) continue end
            r = _gnf_index_Exyz(jE,jx,1,1,Nmx,1,NmE,isFC)
            if scalar
                push!(row,r); push!(coef, C[n] * s^(n-1) * ω[n+1])
            else
                for i in 1:ncross
                    push!(row,r); push!(coef, C[n] * s^(n-1) * ω[n+1,t,i])
                end
            end
        end
        push!(off, Int32(length(row)+1))
    end

    return GNFastClosure(ω0,face,off,row,coef)
end

"""
    sn_fast_applicable(Ndims,Δs,Nmat,𝒪,isFC,Nd,isAdapt;max_cache_bytes)

Whether the optimized SN chain can serve this problem, and why not when it cannot.

Two conditions. The cell-system cache must fit: it holds one factorization per (material,
mesh-width combination, direction), which is bounded only if the mesh is uniform within each
region. And the scheme must not be adaptive — AWD recomputes the closure weights per cell from
the local flux, so the matrix genuinely depends on the voxel and no cache is valid. The
adaptive case falls back to the reference chain.
"""
function sn_fast_applicable(Ndims::Int64,Δs::Vector{Vector{Float64}},Nmat::Int64,
                            𝒪::Vector{Int64},isFC::Bool,Nd::Int64,isAdapt::Bool;
                            max_cache_bytes::Int64=1<<29)

    if Ndims ∉ (1,2,3)
        return false, "the optimized chain covers 1D, 2D and 3D"
    end
    if isAdapt
        return false, "the adaptive (AWD) scheme recomputes the closure weights per cell, " *
                      "which makes the cell system voxel-dependent"
    end

    Nkx = length(unique(Δs[1]))
    Nky = Ndims ≥ 2 ? length(unique(Δs[2])) : 1
    Nkz = Ndims ≥ 3 ? length(unique(Δs[3])) : 1
    Nconf = Nmat*Nkx*Nky*Nkz
    Nm    = isFC ? prod(𝒪) : (1+sum(𝒪.-1))
    bytes = Nd*Nconf*Nm*Nm*8

    if bytes > max_cache_bytes
        return false, "the cell-system cache would need " *
                      "$(round(bytes/2^30,digits=2)) GiB " *
                      "($Nconf material/mesh-width configurations × $Nd directions × " *
                      "$(Nm)² ) — use a mesh that is uniform within each region"
    end
    return true, ""
end

"""
    sn_fast_context(Ndims,Δs,Nmat,Σt,S⁻,S⁺,S,ΔE,𝒪,isFC,isCSD,C,ω,𝒲,Ω,Nd)

Build the optimized SN chain's context for one energy group.

The matrices come from the `sn_*_matrix!` routines, which carry the reference assembly loops
verbatim, and are factorized with `lu!` — once per (material, mesh-width combination,
direction) instead of once per voxel per direction per pass, where the reference pays a full
`𝒮\\Q` including the assembly, the factorization and three allocations.

The six reference kernels do not share one closure shape. `flux_3D_BTE`, `flux_3D_BFP` and
`flux_2D_BFP` index the weights by the *column's* transverse moments only, and reduce to a
common form once the weight arrays are given a uniform `(axis+1, other, other, energy)` layout.
`flux_2D_BTE` and `flux_1D_BFP` instead carry the row's transverse index too and sum the
closure over it — the full transverse coupling — and `flux_1D_BTE` has scalar weights. Hence
the branches below and the two dedicated closure builders.

`Ω` holds the direction-cosine vectors, `Nd` the number of ordinates. Requires `~isAdapt`; see
`sn_fast_applicable`.
"""
function sn_fast_context(Ndims::Int64,Δs::Vector{Vector{Float64}},Nmat::Int64,
                         Σt::Vector{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},
                         S::Array{Float64},ΔE::Float64,𝒪::Vector{Int64},isFC::Bool,
                         isCSD::Bool,C::Vector{Float64},ω::Vector{Array{Float64}},
                         𝒲::Array{Float64},Ω::Vector{Vector{Float64}},Nd::Int64)

    # A genuine 2D (resp. 1D) problem has no z (resp. no y,z) streaming term, so the unused
    # axes are collapsed to a single moment — the 3D tables then reduce exactly to the
    # reference 2D/1D index maps.
    𝒪x = 𝒪[1]
    𝒪y = Ndims ≥ 2 ? 𝒪[2] : 1
    𝒪z = Ndims ≥ 3 ? 𝒪[3] : 1
    𝒪E = isCSD ? 𝒪[4] : 1
    mom = gn_fast_moments(𝒪x,𝒪y,𝒪z,𝒪E,isFC)
    Nm  = mom.Nmom

    Δxv,kx = _gn_fast_distinct(Δs[1])
    Δyv,ky = Ndims ≥ 2 ? _gn_fast_distinct(Δs[2]) : ([1.0],ones(Int32,1))
    Δzv,kz = Ndims ≥ 3 ? _gn_fast_distinct(Δs[3]) : ([1.0],ones(Int32,1))
    Nkx,Nky,Nkz = length(Δxv),length(Δyv),length(Δzv)
    Nconf = Nmat*Nkx*Nky*Nkz
    cells = GNFastCells(Nconf,Nmat,Nkx,Nky,Nkz,kx,ky,kz,Δxv,Δyv,Δzv)

    # `scheme_weights` returns one entry per existing axis, and puts the *energy* first when
    # the CSD term is present: `[ωE, ωx, …]` against `[ωx, …]`. Same convention as the
    # reference kernel calls in the sweeps.
    a0 = isCSD ? 1 : 0
    ωx = ω[a0+1]
    ωy = Ndims ≥ 2 ? ω[a0+2] : zeros(1,1,1)
    ωz = Ndims ≥ 3 ? ω[a0+3] : zeros(1,1,1)
    ωE = isCSD ? ω[1] : zeros(1,1,1,1)

    # Uniform (axis+1, other, other, energy) layout for the kernels whose weights depend on
    # the column's transverse moments only. Unused for 2D BTE and 1D, which are handled apart.
    generic = (Ndims == 3) || (Ndims == 2 && isCSD)
    ωx4 = zeros(1,1,1,1); ωy4 = zeros(1,1,1,1); ωz4 = zeros(1,1,1,1); ωE4 = zeros(1,1,1,1)
    if Ndims == 3 && isCSD                       # already in that layout
        ωx4,ωy4,ωz4,ωE4 = ωx,ωy,ωz,ωE
    elseif Ndims == 3                            # BTE: no energy index
        ωx4 = reshape(ωx, 𝒪x+1, 𝒪y, 𝒪z, 1)
        ωy4 = reshape(ωy, 𝒪y+1, 𝒪x, 𝒪z, 1)
        ωz4 = reshape(ωz, 𝒪z+1, 𝒪x, 𝒪y, 1)
    elseif Ndims == 2 && isCSD                   # no z index
        ωx4 = reshape(ωx, 𝒪x+1, 𝒪y, 1, 𝒪E)
        ωy4 = reshape(ωy, 𝒪y+1, 𝒪x, 1, 𝒪E)
        ωE4 = reshape(ωE, 𝒪E+1, 𝒪x, 𝒪y, 1)
    end

    dirs = Vector{SNFastDir{Nm}}(undef,Nd)
    𝒮 = zeros(Nm,Nm)
    empty_cl = GNFastClosure(Float64[],Int32[],Int32[1],Int32[],Float64[])

    for n in 1:Nd
        μ = Ω[1][n]
        η = Ndims ≥ 2 ? Ω[2][n] : 0.0
        ξ = Ndims ≥ 3 ? Ω[3][n] : 0.0
        isTangential = Ndims == 1 && abs(μ) ≤ 1e-10
        sx = isTangential ? 1.0 : sign(μ)
        sy, sz = sign(η), sign(ξ)
        μk = isTangential ? 0.0 : μ

        LU   = Array{Float64,3}(undef,Nm,Nm,Nconf)
        ipiv = Matrix{Int64}(undef,Nm,Nconf)
        for m in 1:Nmat, ikx in 1:Nkx, iky in 1:Nky, ikz in 1:Nkz
            conf = m + Nmat*((ikx-1) + Nkx*((iky-1) + Nky*(ikz-1)))
            if Ndims == 3
                if isCSD
                    sn_3D_BFP_matrix!(𝒮,μ,η,ξ,Σt[m],S⁺[m],S[m,:],ΔE,Δxv[ikx],Δyv[iky],Δzv[ikz],
                                      𝒪E,𝒪x,𝒪y,𝒪z,C,ωE,ωx,ωy,ωz,𝒲,isFC)
                else
                    sn_3D_BTE_matrix!(𝒮,μ,η,ξ,Σt[m],Δxv[ikx],Δyv[iky],Δzv[ikz],
                                      𝒪x,𝒪y,𝒪z,C,ωx,ωy,ωz,isFC)
                end
            elseif Ndims == 2
                if isCSD
                    sn_2D_BFP_matrix!(𝒮,μ,η,Σt[m],S⁺[m],S[m,:],ΔE,Δxv[ikx],Δyv[iky],
                                      𝒪E,𝒪x,𝒪y,C,ωE,ωx,ωy,𝒲,isFC)
                else
                    sn_2D_BTE_matrix!(𝒮,μ,η,Σt[m],Δxv[ikx],Δyv[iky],𝒪x,𝒪y,C,ωx,ωy,isFC)
                end
            else
                if isCSD
                    sn_1D_BFP_matrix!(𝒮,μ,Σt[m],S⁺[m],S[m,:],ΔE,Δxv[ikx],𝒪E,𝒪x,C,ωE,ωx,𝒲,isFC)
                else
                    sn_1D_BTE_matrix!(𝒮,μ,Σt[m],Δxv[ikx],𝒪x,C,ωx,isFC)
                end
            end
            F = lu!(𝒮)
            LU[:,:,conf]  = F.factors
            ipiv[:,conf] .= F.ipiv
        end

        # Right-hand-side coefficients, transcribed from the reference source loops. The sign
        # is folded in: the reference subtracts these contributions.
        qcx = zeros(Nm,Nkx); qcy = zeros(Nm,Nky); qcz = zeros(Nm,Nkz); qcE = zeros(Nm,Nmat)
        for k in 1:mom.Nact
            jx,jy,jz,jE = mom.ax[k],mom.ay[k],mom.az[k],mom.aE[k]
            wx = generic ? ωx4[1,jy,jz,jE] : (Ndims == 2 ? ωx[1,jy,jy] :
                                              isCSD ? ωx[1,jE,jE] : ωx[1])
            for ikx in 1:Nkx
                hx = abs(μk)/Δxv[ikx]
                qcx[k,ikx] = -C[jx]*hx*(sx^(jx-1)*wx - (-sx)^(jx-1))
            end
            if Ndims ≥ 2
                wy = generic ? ωy4[1,jx,jz,jE] : ωy[1,jx,jx]
                for iky in 1:Nky
                    hy = abs(η)/Δyv[iky]
                    qcy[k,iky] = -C[jy]*hy*(sy^(jy-1)*wy - (-sy)^(jy-1))
                end
            end
            if Ndims ≥ 3
                for ikz in 1:Nkz
                    hz = abs(ξ)/Δzv[ikz]
                    qcz[k,ikz] = -C[jz]*hz*(sz^(jz-1)*ωz4[1,jx,jy,jE] - (-sz)^(jz-1))
                end
            end
            if isCSD
                wE = Ndims == 1 ? ωE[1,jx,jx] : ωE4[1,jx,jy,jz]
                for m in 1:Nmat
                    qcE[k,m] = -C[jE]*((-1)^(jE-1)*S⁺[m]*wE - S⁻[m])
                end
            end
        end

        if generic
            clx = _sn_fast_closure(:x,mom,isFC,ωx4,C,sx)
            cly = Ndims ≥ 2 ? _sn_fast_closure(:y,mom,isFC,ωy4,C,sy) : empty_cl
            clz = Ndims ≥ 3 ? _sn_fast_closure(:z,mom,isFC,ωz4,C,sz) : empty_cl
            clE = isCSD      ? _sn_fast_closure(:E,mom,isFC,ωE4,C,-1.0) : empty_cl
        elseif Ndims == 2
            clx = _sn_closure_2D_bte(:x,mom,isFC,ωx,C,sx)
            cly = _sn_closure_2D_bte(:y,mom,isFC,ωy,C,sy)
            clz = empty_cl; clE = empty_cl
        else
            clx = _sn_closure_1D(:x,mom,isFC,ωx,C,sx,isTangential)
            cly = empty_cl; clz = empty_cl
            clE = isCSD ? _sn_closure_1D(:E,mom,isFC,ωE,C,-1.0,false) : empty_cl
        end

        dirs[n] = SNFastDir{Nm}(Nm,μ ≥ 0 ? 1 : -1,η ≥ 0 ? 1 : -1,ξ ≥ 0 ? 1 : -1,
                                LU,ipiv,qcx,qcy,qcz,qcE,clx,cly,clz,clE)
    end

    return SNFastContext{Nm}(mom,cells,dirs)
end

"""
    SNFastScratch

Per-group scratch for the optimized SN chain: everything the reference reallocates on every
pass, or slices out of a larger array on every direction.

The sweeps of a block of `Nblk` directions run concurrently, so each slot of the block owns a
workspace, its incoming (`src*`) and outgoing (`out*`) face buffers, and a stripe of `𝚿` — the
in-cell moments of that direction, laid out flat so the block's discrete-to-moment transform is
a single `gemm` instead of `Np × Nmom` scattered updates in every cell of every direction.

`Mx`/`Dx` (and their y/z counterparts) hold the half-range transforms of every direction, built
once from `Mn_surf`/`Dn_surf`/`n_to_n⁺` instead of being sliced out per direction and per pass;
`sfc` holds the surface-source intensities as concrete arrays, where the reference carries them
in an abstractly-typed `Union{Float64,Array{Float64}}` container that defeats inference at
every access.
"""
struct SNFastScratch
    Nblk ::Int64
    M    ::Int64                          # Nmom*Nvox — stride between block slots of 𝚿
    ws   ::Vector{GNFastWorkspace}
    𝚿    ::Vector{Float64}                # (M, Nblk)
    srcx ::Vector{Vector{Float64}}
    srcy ::Vector{Vector{Float64}}
    srcz ::Vector{Vector{Float64}}
    outx ::Vector{Vector{Float64}}
    outy ::Vector{Vector{Float64}}
    outz ::Vector{Vector{Float64}}
    Mnn  ::Vector{Vector{Float64}}        # Mn[n,:] of the slot's direction
    Mx   ::Matrix{Float64}                # (Np_surf,Nd) incoming half-range transform
    My   ::Matrix{Float64}
    Mz   ::Matrix{Float64}
    Dx   ::Matrix{Float64}                # (Np_surf,Nd) outgoing half-range transform
    Dy   ::Matrix{Float64}
    Dz   ::Matrix{Float64}
    sfc  ::Vector{Array{Float64,3}}       # (Np_source,·,·) per geometry surface
    Nmf  ::Vector{Int64}                  # face moment counts per axis
    Nfc  ::Vector{Int64}                  # face cell counts per axis
end

"""
    sn_fast_scratch(Ndims,Ns,Nmom,Nm,Nd,Np,Np_surf,Np_source,Mn,Mn_surf,Dn_surf,n_to_n⁺,
    sources;max_block_bytes)

Allocate the per-group scratch of the optimized SN chain and precompute the per-direction
half-range transforms.

`Nblk` is the number of directions swept per batch: their discrete-to-moment transform is
done as a single `gemm` over the batch rather than one `gemv` per direction. It is capped by
memory alone — the block holds `Nblk × Nmom × Nvox` in-cell moments, so on a large mesh it is
allowed to grow only to a fraction of what the flux array itself occupies.
"""
function sn_fast_scratch(Ndims::Int64,Ns::Vector{Int64},Nmom::Int64,Nm::Vector{Int64},
                         Nd::Int64,Np::Int64,Np_surf::Int64,Np_source::Int64,
                         Mn::Array{Float64,2},Mn_surf::Vector{Array{Float64}},
                         Dn_surf::Vector{Array{Float64}},n_to_n⁺::Vector{Vector{Int64}},
                         sources::Array{Union{Array{Float64},Float64}};
                         max_block_bytes::Int64=1<<29)

    Nx,Ny,Nz = Ns[1],Ns[2],Ns[3]
    Nvox = Nx*Ny*Nz
    M    = Nmom*Nvox
    Nblk = min(Nd, max(1, div(max_block_bytes, max(M,1)*8)))

    Nmf = [Nm[1], Ndims ≥ 2 ? Nm[2] : 1, Ndims ≥ 3 ? Nm[3] : 1]
    Nfc = [Ny*Nz, Nx*Nz, Nx*Ny]

    ws   = [GNFastWorkspace(Ndims,Nmom,1,Nmf,Ns) for _ in 1:Nblk]
    srcx = [zeros(Nmf[1]*Nfc[1]) for _ in 1:Nblk]
    srcy = [zeros(Ndims ≥ 2 ? Nmf[2]*Nfc[2] : 0) for _ in 1:Nblk]
    srcz = [zeros(Ndims ≥ 3 ? Nmf[3]*Nfc[3] : 0) for _ in 1:Nblk]
    outx = [zeros(Nmf[1]*Nfc[1]) for _ in 1:Nblk]
    outy = [zeros(Ndims ≥ 2 ? Nmf[2]*Nfc[2] : 0) for _ in 1:Nblk]
    outz = [zeros(Ndims ≥ 3 ? Nmf[3]*Nfc[3] : 0) for _ in 1:Nblk]
    Mnn  = [zeros(Np) for _ in 1:Nblk]

    # Half-range transforms of every direction, transcribed from the reference's per-direction
    # selection: the incoming face of an axis is whichever of its two half ranges the direction
    # belongs to, and a direction tangent to the axis gets none.
    function halfrange(a::Int64)
        Mh = zeros(Np_surf,Nd); Dh = zeros(Np_surf,Nd)
        if Ndims < a return Mh,Dh end
        for n in 1:Nd
            n⁻ = n_to_n⁺[2a-1][n]; n⁺ = n_to_n⁺[2a][n]
            if n⁻ != 0
                Mh[:,n] .= Mn_surf[2a-1][n⁻,:]; Dh[:,n] .= Dn_surf[2a-1][:,n⁻]
            elseif n⁺ != 0
                Mh[:,n] .= Mn_surf[2a][n⁺,:];   Dh[:,n] .= Dn_surf[2a][:,n⁺]
            end
        end
        return Mh,Dh
    end
    Mx,Dx = halfrange(1); My,Dy = halfrange(2); Mz,Dz = halfrange(3)

    # Surface sources as concrete arrays
    dims = [(Ny,Nz),(Ny,Nz),(Nx,Nz),(Nx,Nz),(Nx,Ny),(Nx,Ny)]
    sfc  = Array{Float64,3}[]
    for f in 1:2*Ndims
        n1,n2 = dims[f]
        a = zeros(Np_source,n1,n2)
        for p in 1:Np_source
            sp = sources[p,f]
            if sp isa Array
                for i2 in 1:n2, i1 in 1:n1; a[p,i1,i2] = sp[i1,i2] end
            else
                fill!(view(a,p,:,:), sp)
            end
        end
        push!(sfc,a)
    end
    for f in 2*Ndims+1:6; push!(sfc, zeros(0,0,0)) end

    return SNFastScratch(Nblk,M,ws,zeros(M*Nblk),srcx,srcy,srcz,outx,outy,outz,Mnn,
                         Mx,My,Mz,Dx,Dy,Dz,sfc,Nmf,Nfc)
end
