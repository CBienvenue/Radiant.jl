"""
    first_collision_source!(source::Source, ss::Surface_Source)

Split the transport of a surface source into its analytic uncollided flux ψᵤ and
the collided flux solved by the transport solver. The first-collision source
`Q_fc = Σ_g Σs_ℓ(g→g')·φᵤ` (in-group term included) is added to
`source.volume_sources`, and the uncollided moments are stored in
`source.uncollided_flux` (and `source.uncollided_flux_cutoff` for CSD solvers) to
be added to the solver flux after the solve. Handles monodirectional
(`set_direction`) and distributed (`set_angular_moments`) sources, and both
uncollided models for BFP (`set_uncollided_model`). See TR-09 for the derivation.

Available in 1D Cartesian geometry with void x-boundaries, for the GN and SN
solvers (BTE/BFP/BCSD) and the X-/X+ faces.

# Input Argument(s)
- `source::Source` : source structure (updated in place).
- `ss::Surface_Source` : surface source with `fcs = true`.

# Output Argument(s)
N/A

# Reference(s)
- Bienvenue (2026), TR-09 (First Collision Source).

"""
function first_collision_source!(source::Source, ss::Surface_Source)

    geometry = source.geometry
    cross_sections = source.cross_sections
    solver = source.solver
    # Resolve the cross-section-library particle by tag (Fixed_Sources deepcopies
    # the user source, so ss.particle is not the library object).
    if get_tag(ss.particle) ∉ get_tag.(cross_sections.particles) error(string("No cross sections available for ",get_type(ss.particle)," particle.")) end
    particle = cross_sections.particles[findfirst(x -> get_tag(x) == get_tag(ss.particle),cross_sections.particles)]

    #----
    # Guards
    #----
    Ndims = geometry.get_dimension()
    if lowercase(geometry.get_type()) != "cartesian" || Ndims ∉ [1,2,3]
        error("The first-collision source treatment is only available in 1D/2D/3D cartesian geometries.")
    end
    if any(x -> x != 0, geometry.get_boundary_conditions())
        error("The first-collision source treatment requires void boundaries (the analytic uncollided flux is neither reflected nor wrapped).")
    end
    solver_type,is_CSD = solver.get_solver_type()
    if solver_type ∉ [1,2,3]
        error("The first-collision source treatment is only available for the BTE, BFP and BCSD solvers.")
    end
    if ismissing(ss.energy_group) error("The energy group of the surface source is not defined.") end
    if ismissing(ss.location) error("The location of the surface source is not defined.") end
    face = uppercase(ss.location)
    valid_faces = ["X-","X+","Y-","Y+","Z-","Z+"][1:2*Ndims]
    if face ∉ valid_faces
        error("The first-collision source face $face is not available in $(Ndims)D.")
    end
    is_delta = ~ismissing(ss.direction)
    if ~is_delta && ismissing(ss.angular_moments)
        error("A first-collision surface source requires either a direction (set_direction) or half-range angular moments (set_angular_moments).")
    end
    if is_delta
        # The beam must be incoming with respect to its face (inward normal · Ω > 0).
        bounds = [geometry.voxels_boundaries[["x","y","z"][a]] for a in 1:Ndims]
        _,n̂,_,_ = fcs_face_geometry(face,Ndims,bounds)
        if sum(n̂ .* ss.direction) < 1e-12
            error("The monodirectional surface source must be incoming (and not grazing) with respect to its face.")
        end
    end

    #----
    # Multi-dimensional dispatch (ray-traced pencil-beam superposition)
    #----
    if Ndims >= 2
        first_collision_source_nd!(source,ss,particle,solver,face,solver_type,is_CSD,Ndims)
        return
    end

    #----
    # Angular basis, geometry and material data
    #----
    Np,pl,pm,_ = fcs_angular_basis(solver,1)
    if Np != size(source.volume_sources,2) error("Inconsistent angular basis size for the first-collision source.") end
    Ng = size(source.volume_sources,1)
    if ss.energy_group > Ng error("The energy group of the surface source exceeds the number of groups.") end
    Nx = geometry.number_of_voxels["x"]
    Δx = geometry.voxels_width["x"]
    mat = geometry.get_material_per_voxel()[:,1,1]
    Σt = cross_sections.get_total(particle)
    g0 = ss.energy_group
    I = ss.intensity

    #----
    # Uncollided flux moments φᵤ[g,p,i] (group-integrated, cell-averaged)
    #----
    φ_u = zeros(Ng,Np,Nx)
    φ_cut = zeros(Np,Nx)
    if ~is_CSD
        if is_delta
            𝒴 = fcs_direction_basis(pl,pm,ss.direction)
            fcs_bte_delta!(φ_u,𝒴,I,abs(ss.direction[1]),face,g0,Σt[g0,:],Δx,mat,Nx,Np)
        else
            fcs_bte_moments!(φ_u,ss.angular_moments,I,face,g0,pl,pm,Σt[g0,:],Δx,mat,Nx,Np)
        end
    else
        Sb = cross_sections.get_boundary_stopping_powers(particle)
        Eb = cross_sections.get_energy_boundaries(particle)
        ΔE = cross_sections.get_energy_width(particle)
        cells = face == "X-" ? collect(1:Nx) : collect(Nx:-1:1)
        # Goudsmit-Saunderson model (default, BFP only): the uncollided columns
        # carry the exact FP angular redistribution, driven by the restricted
        # momentum transfer; the "straight" model (and BCSD, which has no FP
        # term) uses Tmt ≡ 0, for which the walkers coincide.
        if solver_type == 2 && lowercase(ss.uncollided_model) == "goudsmit-saunderson"
            Tmt = cross_sections.get_momentum_transfer(particle)
        else
            Tmt = zeros(Ng,size(Σt,2))
        end
        is_GS = any(x->x!=0,Tmt)
        # Flat emission spectrum within the source group (the solver's own
        # boundary-source convention), resolved by Gauss quadrature in energy.
        N_E = 8
        x_E,w_E = gauss_legendre(N_E)
        E0s = [0.5*(Eb[g0+1]+Eb[g0]) + 0.5*(Eb[g0]-Eb[g0+1])*x_E[kE] for kE in 1:N_E]
        w_E0 = 0.5 .* w_E
        walk! = (𝒴,wI,μ̂,E0) -> is_GS ?
            fcs_csd_column_gs!(φ_u,φ_cut,𝒴,pl,wI,μ̂,E0,g0,cells,Δx,mat,Σt,Tmt,Sb,Eb,ΔE,Np) :
            fcs_csd_column!(φ_u,φ_cut,𝒴,wI,μ̂,E0,g0,cells,Δx,mat,Σt,Sb,Eb,ΔE,Np)
        if is_delta
            𝒴 = fcs_direction_basis(pl,pm,ss.direction)
            for kE in 1:N_E
                walk!(𝒴,I*w_E0[kE],abs(ss.direction[1]),E0s[kE])
            end
        else
            g_mom = ss.angular_moments
            L_src = length(g_mom)-1
            N_μ = 64
            x_gl,w_gl = gauss_legendre(N_μ)
            for k in 1:N_μ
                μ̂ = 0.5*(x_gl[k]+1)
                Rp = half_range_legendre_polynomials_up_to_L(L_src,μ̂)
                ψk = 0.0
                for p′ in 0:L_src ψk += g_mom[p′+1]*Rp[p′+1] end
                wI = I*0.5*w_gl[k]*ψk
                if wI == 0.0 continue end
                𝒴 = fcs_column_basis(pl,pm,face == "X-" ? μ̂ : -μ̂)
                for kE in 1:N_E
                    walk!(𝒴,wI*w_E0[kE],μ̂,E0s[kE])
                end
            end
        end
    end

    #----
    # First-collision source Q_fc[g',p,i] = Σ_g Σs_ℓ(g→g')·φᵤ[g,p,i] (in-group included,
    # scattering_source moment convention) and storage of the uncollided moments
    #----
    Lmax = maximum(pl)
    Σs = cross_sections.get_scattering(particle,particle,Lmax)
    for i in range(1,Nx)
        m = mat[i]
        for p in range(1,Np)
            l1 = pl[p]+1
            for g in range(1,Ng)
                if φ_u[g,p,i] == 0.0 continue end
                for g′ in range(1,Ng)
                    Σs_l = Σs[m,g,g′,l1]
                    if Σs_l != 0.0 source.volume_sources[g′,p,1,i,1,1] += Σs_l*φ_u[g,p,i] end
                end
            end
        end
    end
    source.uncollided_flux[:,:,1,:,1,1] .+= φ_u
    source.uncollided_flux_cutoff[:,1,:,1,1] .+= φ_cut

    # Source normalization factor (the 1D face has unit measure; the distributed
    # zeroth moment g₀ is the total incident angular flux).
    ss.normalization_factor += is_delta ? I : I*ss.angular_moments[1]

end

"""
    first_collision_source_nd!(source,ss,particle,solver,face,solver_type,is_CSD,Ndims)

Multi-dimensional (2D/3D) first-collision source, as a ray-traced superposition of
pencil beams: the source face patch is subdivided into launch points, the incident
distribution into directions (a single one for a δ beam; an incoming-hemisphere
quadrature for a distributed source), and each (launch point, direction) ray is
deposited by `fcs_ray_bte!` (BTE) or `fcs_ray_csd!` (BFP/BCSD, Goudsmit--Saunderson
or straight). The uncollided moments and the first-collision source are stored
exactly as in 1D. See TR-09.

# Input Argument(s)
- `source::Source`, `ss::Surface_Source`, `particle::Particle`, `solver::Solver`.
- `face::String`, `solver_type::Int64`, `is_CSD::Bool`, `Ndims::Int64`.

# Output Argument(s)
N/A

"""
function first_collision_source_nd!(source::Source,ss::Surface_Source,particle,solver,face::String,solver_type::Int64,is_CSD::Bool,Ndims::Int64)

    geometry = source.geometry
    cross_sections = source.cross_sections
    axname = ["x","y","z"]
    Np,pl,pm,_ = fcs_angular_basis(solver,Ndims)
    if Np != size(source.volume_sources,2) error("Inconsistent angular basis size for the first-collision source.") end
    Ng = size(source.volume_sources,1)
    g0 = ss.energy_group
    if g0 > Ng error("The energy group of the surface source exceeds the number of groups.") end
    Ns = [geometry.number_of_voxels[axname[a]] for a in 1:Ndims]
    Nx = Ns[1]; Ny = Ndims ≥ 2 ? Ns[2] : 1; Nz = Ndims ≥ 3 ? Ns[3] : 1
    bounds = [geometry.voxels_boundaries[axname[a]] for a in 1:Ndims]
    widths = [geometry.voxels_width[axname[a]] for a in 1:Ndims]
    mat3 = geometry.get_material_per_voxel()
    Σt = cross_sections.get_total(particle)
    R = ss.get_ray_refinement()
    Nμ,Nφ = ss.get_uncollided_angular_order()

    _,n̂,_,_ = fcs_face_geometry(face,Ndims,bounds)
    face_samples = fcs_face_samples(ss,geometry,face,R,Ndims,bounds)
    ang_samples = fcs_angular_samples(ss,n̂,Nμ,Nφ)

    φ_u = zeros(Ng,Np,Nx,Ny,Nz)
    φ_cut = zeros(Np,Nx,Ny,Nz)

    if ~is_CSD
        for (Ω,wI) in ang_samples
            μ̂ = abs(sum(Ω .* n̂)); 𝒴 = fcs_direction_basis(pl,pm,Ω)
            for (r0,ΔA) in face_samples
                segs,sin_p = fcs_ray_trace(r0,Ω,bounds,Ndims)
                fcs_ray_bte!(φ_u,segs,𝒴,wI*μ̂*ΔA,sin_p,Σt[g0,:],widths,mat3,g0,Np,Ndims)
            end
        end
    else
        Sb = cross_sections.get_boundary_stopping_powers(particle)
        Eb = cross_sections.get_energy_boundaries(particle)
        ΔE = cross_sections.get_energy_width(particle)
        if solver_type == 2 && lowercase(ss.uncollided_model) == "goudsmit-saunderson"
            Tmt = cross_sections.get_momentum_transfer(particle)
        else
            Tmt = zeros(Ng,size(Σt,2))
        end
        # Flat emission spectrum within the source group (solver convention).
        N_E = 8
        x_E,w_E = gauss_legendre(N_E)
        E0s = [0.5*(Eb[g0+1]+Eb[g0]) + 0.5*(Eb[g0]-Eb[g0+1])*x_E[kE] for kE in 1:N_E]
        w_E0 = 0.5 .* w_E
        for (Ω,wI) in ang_samples
            μ̂ = abs(sum(Ω .* n̂)); 𝒴 = fcs_direction_basis(pl,pm,Ω)
            for (r0,ΔA) in face_samples
                segs,sin_p = fcs_ray_trace(r0,Ω,bounds,Ndims)
                for kE in 1:N_E
                    fcs_ray_csd!(φ_u,φ_cut,segs,𝒴,pl,wI*w_E0[kE]*μ̂*ΔA,sin_p,E0s[kE],g0,Σt,Tmt,Sb,Eb,ΔE,widths,mat3,Np,Ndims,r0,Ω,bounds)
                end
            end
        end
    end

    # First-collision source Q_fc and storage (in-group included; scattering_source convention).
    Lmax = maximum(pl)
    Σs = cross_sections.get_scattering(particle,particle,Lmax)
    for iz in range(1,Nz), iy in range(1,Ny), ix in range(1,Nx)
        m = mat3[ix,iy,iz]
        for p in range(1,Np)
            l1 = pl[p]+1
            for g in range(1,Ng)
                if φ_u[g,p,ix,iy,iz] == 0.0 continue end
                for g′ in range(1,Ng)
                    Σs_l = Σs[m,g,g′,l1]
                    if Σs_l != 0.0 source.volume_sources[g′,p,1,ix,iy,iz] += Σs_l*φ_u[g,p,ix,iy,iz] end
                end
            end
        end
    end
    source.uncollided_flux[:,:,1,:,:,:] .+= φ_u
    source.uncollided_flux_cutoff[:,1,:,:,:] .+= φ_cut

    # Source normalization factor: intensity × illuminated face area (× g₀ for a
    # distributed source), matching the boundary-treatment convention.
    area = sum(ΔA for (_,ΔA) in face_samples)
    g0factor = ismissing(ss.angular_moments) ? 1.0 : ss.angular_moments[1]
    ss.normalization_factor += ss.intensity*area*g0factor

end

"""
    fcs_angular_basis(solver::Solver,Ndims::Int64)

Resolve the solver's full-range angular moment basis for the first-collision
source: number of moments `Np`, Legendre orders `pl`, azimuthal orders `pm` and
the basis type.

# Input Argument(s)
- `solver::Solver` : transport solver (GN or SN).
- `Ndims::Int64` : geometry dimension.

# Output Argument(s)
- `Np::Int64` : number of angular moments.
- `pl::Vector{Int64}` : Legendre order per moment.
- `pm::Vector{Int64}` : azimuthal order per moment.
- `is_SPH::Bool` : spherical-harmonics (true) or Legendre-type (false) basis.

"""
function fcs_angular_basis(solver::Solver,Ndims::Int64)
    if solver isa GN
        L = solver.get_legendre_order()
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        if polynomial_basis == "legendre"
            Np = L+1
            pl = collect(0:L)
            pm = zeros(Int64,Np)
            is_SPH = false
        else
            Np = spherical_harmonics_number_basis(L)
            pl,pm = spherical_harmonics_indices(L)
            is_SPH = true
        end
    elseif solver isa SN
        Qdims = solver.get_quadrature_dimension(Ndims)
        Ω,w = quadrature(solver.quadrature_order,solver.quadrature_type,Qdims)
        if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
        Np,_,_,pl,pm = angular_polynomial_basis(Ω,w,solver.get_legendre_order(),solver.get_angular_boltzmann(),Qdims)
        is_SPH = false
    else
        error("The first-collision source treatment is only available with the GN and SN solvers.")
    end
    return Np,pl,pm,is_SPH
end

"""
    fcs_direction_basis(pl::Vector{Int64},pm::Vector{Int64},Ω::Vector{Float64})

Full-range moment values `𝒴_p` of the angular δ-function at direction `Ω`:
`R_ℓm(μₛ,ϕₛ)` whenever the basis carries azimuthal orders (GN spherical harmonics,
SN 2D/3D), reducing to `P_ℓ(μₛ)` for the pure-Legendre 1D bases (pm ≡ 0).

# Input Argument(s)
- `pl::Vector{Int64}`, `pm::Vector{Int64}` : moment orders.
- `Ω::Vector{Float64}` : direction cosines [μ,η,ξ].

# Output Argument(s)
- `𝒴::Vector{Float64}` : basis values per moment.

"""
function fcs_direction_basis(pl::Vector{Int64},pm::Vector{Int64},Ω::Vector{Float64})
    Np = length(pl)
    Lmax = maximum(pl)
    𝒴 = zeros(Np)
    if any(m -> m != 0, pm)
        ϕs = atan(Ω[3],Ω[2])
        Rlm = real_spherical_harmonics_up_to_L(Lmax,Ω[1],ϕs)
        for p in range(1,Np)
            𝒴[p] = Rlm[pl[p]+1][pl[p]+pm[p]+1]
        end
    else
        Pls = legendre_polynomials_up_to_L(Lmax,Ω[1])
        for p in range(1,Np)
            𝒴[p] = Pls[pl[p]+1]
        end
    end
    return 𝒴
end

"""
    fcs_column_basis(pl::Vector{Int64},pm::Vector{Int64},μ::Float64)

Full-range moment values `P_ℓ(μ)` of an azimuthally-symmetric column at direction
cosine `μ`; only the m = 0 slots are filled.

# Input Argument(s)
- `pl::Vector{Int64}`, `pm::Vector{Int64}` : moment orders.
- `μ::Float64` : direction cosine.

# Output Argument(s)
- `𝒴::Vector{Float64}` : basis values per moment.

"""
function fcs_column_basis(pl::Vector{Int64},pm::Vector{Int64},μ::Float64)
    Np = length(pl)
    𝒴 = zeros(Np)
    Pls = legendre_polynomials_up_to_L(maximum(pl),μ)
    for p in range(1,Np)
        if pm[p] == 0
            𝒴[p] = Pls[pl[p]+1]
        end
    end
    return 𝒴
end

"""
    fcs_bte_delta!(φ_u,𝒴,I,μ̂,face,g0,Σt_g,Δx,mat,Nx,Np)

Closed-form cell-averaged uncollided flux moments of a monodirectional boundary
source in a BTE solve, `φᵤ_p(i) = I·𝒴_p·Wᵢ`, where `Wᵢ` is the cell-averaged
attenuation `e^(-τ/μ̂)` and τ the optical depth from the source face (TR-09).

# Input Argument(s)
- `φ_u::Array{Float64,3}` : uncollided flux moments `[Ng,Np,Nx]` (updated in place).
- `𝒴::Vector{Float64}` : moment-basis values at the source direction.
- `I::Float64` : source intensity.
- `μ̂::Float64` : |direction cosine| of the source.
- `face::String` : source face ("X-" or "X+").
- `g0::Int64` : emission group.
- `Σt_g::Vector{Float64}` : total cross section per material in the emission group.
- `Δx::Vector{Float64}`, `mat::Vector{Int64}`, `Nx::Int64`, `Np::Int64` : mesh data.

# Output Argument(s)
N/A

"""
function fcs_bte_delta!(φ_u::Array{Float64,3},𝒴::Vector{Float64},I::Float64,μ̂::Float64,face::String,g0::Int64,Σt_g::Vector{Float64},Δx::Vector{Float64},mat::Vector{Int64},Nx::Int64,Np::Int64)
    cells = face == "X-" ? (1:Nx) : (Nx:-1:1)
    τ = 0.0
    for i in cells
        Σ = Σt_g[mat[i]]
        W = Σ > 0 ? μ̂/(Σ*Δx[i])*(exp(-τ/μ̂) - exp(-(τ+Σ*Δx[i])/μ̂)) : exp(-τ/μ̂)
        for p in range(1,Np)
            if 𝒴[p] != 0.0 φ_u[g0,p,i] += I*𝒴[p]*W end
        end
        τ += Σ*Δx[i]
    end
end

"""
    fcs_bte_moments!(φ_u,g_mom,I,face,g0,pl,pm,Σt_g,Δx,mat,Nx,Np)

Exact cell-averaged uncollided flux moments of a distributed boundary source
(half-range moments `g_p`) in a BTE solve. The μ-integrals reduce to generalized
exponential integrals `E_n` in closed form via `cp_ray_reduction` (called with
`w = 1`), promoting to BigFloat when the combined polynomial order exceeds 12.
See TR-09.

# Input Argument(s)
- `φ_u::Array{Float64,3}` : uncollided flux moments `[Ng,Np,Nx]` (updated in place).
- `g_mom::Vector{Float64}` : half-range moments of the incident angular flux.
- `I::Float64` : source intensity.
- `face::String` : source face ("X-" or "X+").
- `g0::Int64` : emission group.
- `pl::Vector{Int64}`, `pm::Vector{Int64}` : moment basis orders.
- `Σt_g::Vector{Float64}` : total cross section per material in the emission group.
- `Δx::Vector{Float64}`, `mat::Vector{Int64}`, `Nx::Int64`, `Np::Int64` : mesh data.

# Output Argument(s)
N/A

"""
function fcs_bte_moments!(φ_u::Array{Float64,3},g_mom::Vector{Float64},I::Float64,face::String,g0::Int64,pl::Vector{Int64},pm::Vector{Int64},Σt_g::Vector{Float64},Δx::Vector{Float64},mat::Vector{Int64},Nx::Int64,Np::Int64)
    Lmax = maximum(pl)
    L_src = length(g_mom)-1
    # Moment slots receiving the order-ℓ value (m = 0 only in the SH basis).
    slots = [Int64[] for l in 0:Lmax]
    for p in range(1,Np)
        if pm[p] == 0 push!(slots[pl[p]+1],p) end
    end
    T = (Lmax + L_src > 12) ? BigFloat : Float64
    compute = function()
        c_leg = [cp_legendre_polynomial_coefficients(l,T) for l in 0:Lmax]
        c_half = [T(sqrt(2*p′+1)) .* cp_shifted_legendre_polynomial_coefficients(p′,T) for p′ in 0:L_src]
        c = [[cp_polynomial_product(c_leg[l+1],c_half[p′+1]) for p′ in 0:L_src] for l in 0:Lmax]
        cells = face == "X-" ? (1:Nx) : (Nx:-1:1)
        τ = zero(T)
        for i in cells
            Σi = T(Σt_g[mat[i]])
            Δxi = T(Δx[i])
            for l in range(0,Lmax)
                if isempty(slots[l+1]) continue end
                val = zero(T)
                for p′ in range(0,L_src)
                    if g_mom[p′+1] != 0.0
                        val += T(g_mom[p′+1])*cp_ray_reduction(c[l+1][p′+1],Σi,Δxi,τ,1)
                    end
                end
                w = I*Float64(val)/Δx[i]
                if face == "X+" && isodd(l) w = -w end
                for p in slots[l+1]
                    φ_u[g0,p,i] += w
                end
            end
            τ += Σi*Δxi
        end
    end
    if T == BigFloat
        setprecision(BigFloat,64+8*(Lmax+L_src)) do
            compute()
        end
    else
        compute()
    end
end

"""
    fcs_csd_column_gs!(φ_u,φ_cut,𝒴,pl,wI,μ̂,E0,g0,cells,Δx,mat,Σt,Tmt,Sb,Eb,ΔE,Np)

Goudsmit--Saunderson uncollided column for BFP solves: the event walker of
`fcs_csd_column!` with the exact Fokker--Planck angular redistribution added along
the column pathlength. The order-ℓ moments damp as `e^(-ℓ(ℓ+1)·τ_T)`,
`τ_T = ∫T ds` (Radiant convention `λ_ℓ = ℓ(ℓ+1)·T`, matching
`fokker_planck_galerkin`), giving the per-segment effective attenuation
`Σ̄ + λ_ℓ·T_gm`; the depth advances at the exact mean cosine
`⟨μ⟩ = μ̂₀·e^(-2τ_T)` (first Lewis moment). Only the longitudinal straggling is
neglected. With `Tmt ≡ 0` it reduces to `fcs_csd_column!`. See TR-09.

# Input Argument(s)
Same as `fcs_csd_column!`, plus:
- `pl::Vector{Int64}` : Legendre order per moment slot.
- `Tmt::Matrix{Float64}` : restricted momentum transfer `[Ng,Nmat]`.

# Output Argument(s)
N/A

"""
function fcs_csd_column_gs!(φ_u::Array{Float64,3},φ_cut::Matrix{Float64},𝒴::Vector{Float64},pl::Vector{Int64},wI::Float64,μ̂::Float64,E0::Float64,g0::Int64,cells::Vector{Int64},Δx::Vector{Float64},mat::Vector{Int64},Σt::Matrix{Float64},Tmt::Matrix{Float64},Sb::Matrix{Float64},Eb::Vector{Float64},ΔE::Vector{Float64},Np::Int64)
    Ng = length(ΔE)
    Lmax = maximum(pl)
    fac = zeros(Lmax+1)                    # per-order segment deposit factors
    g = g0
    E = E0
    T_att = 1.0
    τ_T = 0.0
    ic = 1
    i = cells[ic]
    d_cell = Δx[i]                         # depth remaining in the current cell
    while g ≤ Ng && T_att > 1e-300
        m = mat[i]
        E_lo = Eb[g+1]
        k = (Sb[g,m] - Sb[g+1,m])/(Eb[g] - E_lo)
        S_E = Sb[g+1,m] + k*(E - E_lo)
        Δs_grp = abs(k)*ΔE[g] > 1e-12*S_E ? log(S_E/Sb[g+1,m])/k : (E - E_lo)/S_E
        T_gm = Tmt[g,m]
        # Slant path to the cell edge through the mean-depth map
        # dx = μ̂₀·D₁(s)·ds = μ̂₀·e^(-2τ_T)ds (λ₁ = 1·2 in the Radiant convention).
        if T_gm > 0
            arg = 1.0 - 2*T_gm*d_cell*exp(2*τ_T)/μ̂
            Δs_cell = arg > 0 ? -log(arg)/(2*T_gm) : Inf
        else
            Δs_cell = d_cell*exp(2*τ_T)/μ̂
        end
        crosses_group = Δs_grp ≤ Δs_cell
        Δs = crosses_group ? Δs_grp : Δs_cell
        Σ̄ = Σt[g,m]
        # Per-order deposit factors with the effective attenuation Σ̄ + λ_ℓ·T_gm.
        for l in range(0,Lmax)
            Σ_eff = Σ̄ + l*(l+1)*T_gm
            fac[l+1] = exp(-l*(l+1)*τ_T) * (Σ_eff*Δs > 1e-12 ? (1-exp(-Σ_eff*Δs))/Σ_eff : Δs)
        end
        w_seg = wI*(μ̂/Δx[i])*T_att
        for p in range(1,Np)
            if 𝒴[p] != 0.0 φ_u[g,p,i] += w_seg*𝒴[p]*fac[pl[p]+1] end
        end
        # Advance the column state.
        T_att *= exp(-Σ̄*Δs)
        Δd = T_gm > 0 ? μ̂*exp(-2*τ_T)*(1-exp(-2*T_gm*Δs))/(2*T_gm) : μ̂*exp(-2*τ_T)*Δs
        τ_T += T_gm*Δs
        if crosses_group
            d_cell -= Δd
            g += 1
            E = E_lo
            if g > Ng
                # The column reaches the cutoff energy inside cell i.
                w_cut = wI*μ̂*T_att/Sb[Ng+1,m]*ΔE[Ng]/Δx[i]
                for p in range(1,Np)
                    if 𝒴[p] != 0.0 φ_cut[p,i] += w_cut*𝒴[p]*exp(-pl[p]*(pl[p]+1)*τ_T) end
                end
                return
            end
        else
            if abs(k)*ΔE[g] > 1e-12*S_E
                E = E_lo + (S_E*exp(-k*Δs) - Sb[g+1,m])/k
            else
                E -= S_E*Δs
            end
            ic += 1
            if ic > length(cells) return end   # transmitted through the far face
            i = cells[ic]
            d_cell = Δx[i]
        end
    end
end

"""
    fcs_csd_column!(φ_u,φ_cut,𝒴,wI,μ̂,E0,g0,cells,Δx,mat,Σt,Sb,Eb,ΔE,Np)

Exact multigroup projection of one straight uncollided column in a CSD solve. The
walker advances event by event (cell and group-boundary crossings), with energy
and pathlength locked by the range relation (stopping power linear in energy
between the group-boundary values `Sb`); it deposits the group-integrated
cell-averaged flux into `φ_u` and, at the cutoff energy, into `φ_cut` in the
solver's `𝚽cutoff` convention. See TR-09.

# Input Argument(s)
- `φ_u::Array{Float64,3}` : uncollided flux moments `[Ng,Np,Nx]` (updated in place).
- `φ_cut::Matrix{Float64}` : cutoff flux moments `[Np,Nx]` (updated in place).
- `𝒴::Vector{Float64}` : moment-basis values at the column direction.
- `wI::Float64` : column weight (intensity, or quadrature-weighted intensity).
- `μ̂::Float64` : |direction cosine| of the column.
- `E0::Float64` : emission energy.
- `g0::Int64` : emission group.
- `cells::Vector{Int64}` : cell indices in traversal order.
- `Δx::Vector{Float64}`, `mat::Vector{Int64}` : mesh data.
- `Σt::Matrix{Float64}` : catastrophic total cross sections `[Ng,Nmat]`.
- `Sb::Matrix{Float64}` : group-boundary stopping powers `[Ng+1,Nmat]`.
- `Eb::Vector{Float64}` : group energy boundaries `[Ng+1]`.
- `ΔE::Vector{Float64}` : group energy widths `[Ng]`.
- `Np::Int64` : number of angular moments.

# Output Argument(s)
N/A

"""
function fcs_csd_column!(φ_u::Array{Float64,3},φ_cut::Matrix{Float64},𝒴::Vector{Float64},wI::Float64,μ̂::Float64,E0::Float64,g0::Int64,cells::Vector{Int64},Δx::Vector{Float64},mat::Vector{Int64},Σt::Matrix{Float64},Sb::Matrix{Float64},Eb::Vector{Float64},ΔE::Vector{Float64},Np::Int64)
    Ng = length(ΔE)
    g = g0
    E = E0
    T = 1.0
    ic = 1
    i = cells[ic]
    s_cell = Δx[i]/μ̂                      # slant path remaining in the current cell
    while g ≤ Ng && T > 1e-300
        m = mat[i]
        E_lo = Eb[g+1]
        k = (Sb[g,m] - Sb[g+1,m])/(Eb[g] - E_lo)    # dS/dE within the group (S linear in E)
        S_E = Sb[g+1,m] + k*(E - E_lo)
        Δs_grp = abs(k)*ΔE[g] > 1e-12*S_E ? log(S_E/Sb[g+1,m])/k : (E - E_lo)/S_E
        crosses_group = Δs_grp ≤ s_cell
        Δs = crosses_group ? Δs_grp : s_cell
        Σ̄ = Σt[g,m]
        att = exp(-Σ̄*Δs)
        seg = Σ̄ > 0 ? (1-att)/Σ̄ : Δs
        w_seg = wI*(μ̂/Δx[i])*T*seg
        for p in range(1,Np)
            if 𝒴[p] != 0.0 φ_u[g,p,i] += w_seg*𝒴[p] end
        end
        T *= att
        if crosses_group
            s_cell -= Δs
            g += 1
            E = E_lo
            if g > Ng
                # The column reaches the cutoff energy inside cell i.
                w_cut = wI*T*μ̂/Sb[Ng+1,m]*ΔE[Ng]/Δx[i]
                for p in range(1,Np)
                    if 𝒴[p] != 0.0 φ_cut[p,i] += w_cut*𝒴[p] end
                end
                return
            end
        else
            # Advance to the next cell, updating the energy after the partial in-group path.
            if abs(k)*ΔE[g] > 1e-12*S_E
                E = E_lo + (S_E*exp(-k*Δs) - Sb[g+1,m])/k
            else
                E -= S_E*Δs
            end
            ic += 1
            if ic > length(cells) return end   # transmitted through the far face
            i = cells[ic]
            s_cell = Δx[i]/μ̂
        end
    end
end

"""
    fcs_ray_trace(r0::Vector{Float64},Ω::Vector{Float64},bounds::Vector{Vector{Float64}},Ndims::Int64)

Amanatides--Woo traversal of a ray through a non-uniform Cartesian voxel mesh. The
ray starts at `r0` (on a source face) and travels along the active-axis direction
`Ω[1:Ndims]`; it returns the ordered list of crossed voxels with their geometric
segment length (in active-axis space), until the ray exits the mesh. The active
direction norm `sin_p = |Ω[1:Ndims]|` converts a geometric segment to true
pathlength (`ds_true = ℓ/sin_p`), accounting for out-of-plane travel in 2D.

# Input Argument(s)
- `r0::Vector{Float64}` : ray origin (length ≥ Ndims).
- `Ω::Vector{Float64}` : direction cosines [μ,η,ξ].
- `bounds::Vector{Vector{Float64}}` : voxel boundary array per active axis.
- `Ndims::Int64` : number of active spatial dimensions.

# Output Argument(s)
- `segments::Vector{Tuple{Int64,Int64,Int64,Float64}}` : `(ix,iy,iz,ℓ_geom)` per voxel.
- `sin_p::Float64` : norm of the active-axis direction.

"""
function fcs_ray_trace(r0::Vector{Float64},Ω::Vector{Float64},bounds::Vector{Vector{Float64}},Ndims::Int64)
    d = Ω[1:Ndims]
    sin_p = sqrt(sum(d.^2))
    d = d ./ sin_p                       # unit direction in active-axis space
    N = [length(bounds[a])-1 for a in 1:Ndims]
    cell = zeros(Int64,Ndims)
    step = zeros(Int64,Ndims)
    tnext = fill(Inf,Ndims)
    for a in 1:Ndims
        b = bounds[a]
        if d[a] > 0
            c = clamp(searchsortedlast(b,r0[a]+1e-12),1,N[a]); step[a] = 1
            tnext[a] = (b[c+1]-r0[a])/d[a]
        elseif d[a] < 0
            c = clamp(searchsortedlast(b,r0[a]-1e-12),1,N[a]); step[a] = -1
            tnext[a] = (b[c]-r0[a])/d[a]
        else
            c = clamp(searchsortedlast(b,r0[a]),1,N[a]); step[a] = 0
        end
        cell[a] = c
    end
    segments = Tuple{Int64,Int64,Int64,Float64}[]
    t = 0.0
    while all(1 ≤ cell[a] ≤ N[a] for a in 1:Ndims)
        amin = 1; tcross = tnext[1]
        for a in 2:Ndims
            if tnext[a] < tcross amin = a; tcross = tnext[a] end
        end
        ℓ = tcross - t
        ix = cell[1]; iy = Ndims ≥ 2 ? cell[2] : 1; iz = Ndims ≥ 3 ? cell[3] : 1
        if ℓ > 1e-14 push!(segments,(ix,iy,iz,ℓ)) end
        t = tcross
        cell[amin] += step[amin]
        b = bounds[amin]; c = cell[amin]
        tnext[amin] = (1 ≤ c ≤ N[amin]) ? (step[amin] > 0 ? (b[c+1]-r0[amin])/d[amin] : (b[c]-r0[amin])/d[amin]) : Inf
    end
    return segments,sin_p
end

"""
    fcs_face_geometry(face::String,Ndims::Int64,bounds::Vector{Vector{Float64}})

Geometry of a source face: the beam axis, the inward unit normal `n̂`, the tangent
axes and the beam-axis coordinate of the face plane.

# Input Argument(s)
- `face::String` : source face ("X-","X+","Y-","Y+","Z-","Z+").
- `Ndims::Int64` : number of active spatial dimensions.
- `bounds::Vector{Vector{Float64}}` : voxel boundary array per active axis.

# Output Argument(s)
- `axis::Int64` : beam axis (1=x,2=y,3=z).
- `n̂::Vector{Float64}` : inward unit normal.
- `tangents::Vector{Int64}` : tangent axes.
- `plane::Float64` : beam-axis coordinate of the face.

"""
function fcs_face_geometry(face::String,Ndims::Int64,bounds::Vector{Vector{Float64}})
    axis = face[1] == 'X' ? 1 : face[1] == 'Y' ? 2 : 3
    plus = face[2] == '+'
    n̂ = zeros(3); n̂[axis] = plus ? -1.0 : 1.0
    plane = plus ? bounds[axis][end] : bounds[axis][1]
    tangents = [a for a in 1:Ndims if a != axis]
    return axis,n̂,tangents,plane
end

"""
    fcs_face_samples(ss::Surface_Source,geometry::Geometry,face::String,R::Int64,Ndims::Int64,bounds)

Launch points and face areas of the illuminated source patch, subdivided into
`R^(Ndims-1)` sub-elements per face cell (a face cell is illuminated when its
tangent center lies within `set_boundaries`, the `surface_source.jl` convention).

# Input Argument(s)
- `ss::Surface_Source`, `geometry::Geometry`, `face::String`, `R::Int64`, `Ndims::Int64`.
- `bounds::Vector{Vector{Float64}}` : voxel boundary array per active axis.

# Output Argument(s)
- `samples::Vector{Tuple{Vector{Float64},Float64}}` : `(r0, ΔA_face)` per sub-element.

"""
function fcs_face_samples(ss::Surface_Source,geometry::Geometry,face::String,R::Int64,Ndims::Int64,bounds::Vector{Vector{Float64}})
    axis,_,tangents,plane = fcs_face_geometry(face,Ndims,bounds)
    axname = ["x","y","z"]
    tsub = [Tuple{Float64,Float64}[] for _ in tangents]
    for (ti,a) in enumerate(tangents)
        pos = geometry.voxels_position[axname[a]]
        wid = geometry.voxels_width[axname[a]]
        bnd = bounds[a]
        lo,hi = haskey(ss.boundaries,axname[a]) ? (ss.boundaries[axname[a]][1],ss.boundaries[axname[a]][2]) : (-Inf,Inf)
        for ic in 1:length(pos)
            if pos[ic] < lo || pos[ic] > hi continue end
            for r in 1:R
                subw = wid[ic]/R
                push!(tsub[ti],(bnd[ic]+(r-0.5)*subw,subw))
            end
        end
    end
    samples = Tuple{Vector{Float64},Float64}[]
    if length(tangents) == 1
        for (c,w) in tsub[1]
            r0 = zeros(3); r0[axis] = plane; r0[tangents[1]] = c
            push!(samples,(r0,w))
        end
    else
        for (c1,w1) in tsub[1], (c2,w2) in tsub[2]
            r0 = zeros(3); r0[axis] = plane; r0[tangents[1]] = c1; r0[tangents[2]] = c2
            push!(samples,(r0,w1*w2))
        end
    end
    return samples
end

"""
    fcs_angular_samples(ss::Surface_Source,n̂::Vector{Float64},Nμ::Int64,Nφ::Int64)

Direction samples of the incident distribution: a single `(Ω_s, I)` for a
monodirectional source, or a Gauss(μ̂)×uniform(azimuth) quadrature of the incoming
hemisphere about `n̂` for a distributed source (`wI_k = w_k·I·Σ_p g_p R̄_p(μ̂_k)`).

# Input Argument(s)
- `ss::Surface_Source`, `n̂::Vector{Float64}`.
- `Nμ::Int64`, `Nφ::Int64` : polar and azimuthal quadrature orders (distributed only).

# Output Argument(s)
- `samples::Vector{Tuple{Vector{Float64},Float64}}` : `(Ω, wI)` per direction.

"""
function fcs_angular_samples(ss::Surface_Source,n̂::Vector{Float64},Nμ::Int64,Nφ::Int64)
    if ismissing(ss.angular_moments)
        return [(ss.direction,ss.intensity)]
    end
    g = ss.angular_moments; Lsrc = length(g)-1
    # Tangent frame (n̂ is along a coordinate axis).
    ax = argmax(abs.(n̂))
    e1 = zeros(3); e1[mod1(ax+1,3)] = 1.0
    e2 = zeros(3); e2[mod1(ax+2,3)] = 1.0
    xμ,wμ = gauss_legendre(Nμ)
    samples = Tuple{Vector{Float64},Float64}[]
    for k in 1:Nμ
        μ̂ = 0.5*(xμ[k]+1); wk = 0.5*wμ[k]
        sinθ = sqrt(max(0.0,1-μ̂^2))
        Pl = legendre_polynomials_up_to_L(Lsrc,2*μ̂-1)
        ψ = 0.0
        for p in 0:Lsrc ψ += g[p+1]*sqrt(2*p+1)*Pl[p+1] end
        # The moments g describe the azimuthally-INTEGRATED incident flux: the
        # per-steradian distribution is ψ_az/2π, so the azimuth is averaged.
        for j in 1:Nφ
            α = 2π*(j-0.5)/Nφ
            Ω = μ̂ .* n̂ .+ (sinθ*cos(α)) .* e1 .+ (sinθ*sin(α)) .* e2
            push!(samples,(Ω,wk*ψ*ss.intensity/Nφ))
        end
    end
    return samples
end

"""
    fcs_ray_bte!(φ_u,segs,𝒴,wA,sin_p,Σt_g,widths,mat3,g0,Np,Ndims)

Deposit the BTE uncollided flux moments of one pencil beam (ray `segs`) into the
`(ix,iy,iz)` voxels: `φᵤ[g0,p,c] += (wA/V_c)·𝒴_p·∫_seg e^(-τ) ds_true`, with true
path `ds = ℓ/sin_p` and `wA = wI·μ̂·ΔA_face`.

# Input Argument(s)
- `φ_u::Array{Float64,5}` : `[Ng,Np,Nx,Ny,Nz]` (updated in place).
- `segs` : traced `(ix,iy,iz,ℓ)` segments. `𝒴::Vector{Float64}` : moment values at Ω.
- `wA::Float64` : ray weight `wI·μ̂·ΔA_face`. `sin_p::Float64` : active-direction norm.
- `Σt_g::Vector{Float64}` : total cross section per material (emission group).
- `widths::Vector{Vector{Float64}}`, `mat3::Array{Int64,3}`, `g0::Int64`, `Np::Int64`, `Ndims::Int64`.

# Output Argument(s)
N/A

"""
function fcs_ray_bte!(φ_u::Array{Float64,5},segs::Vector{Tuple{Int64,Int64,Int64,Float64}},𝒴::Vector{Float64},wA::Float64,sin_p::Float64,Σt_g::Vector{Float64},widths::Vector{Vector{Float64}},mat3::Array{Int64,3},g0::Int64,Np::Int64,Ndims::Int64)
    τ = 0.0
    for (ix,iy,iz,ℓ) in segs
        m = mat3[ix,iy,iz]; Σ = Σt_g[m]; ds = ℓ/sin_p
        τ_out = τ + Σ*ds
        seg = Σ > 0 ? (exp(-τ)-exp(-τ_out))/Σ : ds*exp(-τ)
        Vc = widths[1][ix]*(Ndims ≥ 2 ? widths[2][iy] : 1.0)*(Ndims ≥ 3 ? widths[3][iz] : 1.0)
        w = wA/Vc*seg
        for p in range(1,Np)
            if 𝒴[p] != 0.0 φ_u[g0,p,ix,iy,iz] += w*𝒴[p] end
        end
        τ = τ_out
    end
end

"""
    fcs_ray_csd!(φ_u,φ_cut,segs,𝒴,pl,wA,sin_p,E0,g0,Σt,Tmt,Sb,Eb,ΔE,widths,mat3,Np,Ndims,r0,Ω,bounds)

Deposit the CSD/BFP uncollided flux of one pencil along the traced ray, the
multi-dimensional generalization of `fcs_csd_column_gs!`: event walker over voxel
crossings (from `segs`) and group crossings (CSD range relation), with the exact
Goudsmit--Saunderson angular damping (`e^(-ℓ(ℓ+1)τ_T)`, effective attenuation
`Σ̄+ℓ(ℓ+1)T`) and the mean geometric advance compressed by `D₁=e^(-2τ_T)`. Each
deposit is spread spatially by a Gaussian kernel carrying the exact Lewis second
moments of the pencil (`fcs_lewis_update`): longitudinal variance `⟨z²⟩-z̄²` along
the beam axis and transverse variance `⟨x²⟩` per perpendicular axis, projected on
the grid axes (`σ_a² = σ⊥²+(σ_z²-σ⊥²)·Ω_a²`, `fcs_splat_windows`). With `Tmt ≡ 0`
(straight, BCSD) the kernel is degenerate and the deposit is single-voxel. Cutoff
deposited in the solver convention, spread by the kernel at the cutoff pathlength.

# Input Argument(s)
- `φ_u::Array{Float64,5}`, `φ_cut::Array{Float64,4}` : updated in place.
- `segs`, `𝒴`, `pl` : ray segments, moment values at Ω, Legendre orders.
- `wA::Float64` (`wI·μ̂·ΔA_face`), `sin_p::Float64`, `E0::Float64`, `g0::Int64`.
- `Σt,Tmt::Matrix{Float64}` (`[Ng,Nmat]`), `Sb::Matrix{Float64}`, `Eb,ΔE::Vector{Float64}`.
- `widths`, `mat3`, `Np`, `Ndims`.
- `r0::Vector{Float64}`, `Ω::Vector{Float64}`, `bounds::Vector{Vector{Float64}}` :
  ray origin, direction and voxel boundaries per active axis (kernel placement).

# Output Argument(s)
N/A

"""
function fcs_ray_csd!(φ_u::Array{Float64,5},φ_cut::Array{Float64,4},segs::Vector{Tuple{Int64,Int64,Int64,Float64}},𝒴::Vector{Float64},pl::Vector{Int64},wA::Float64,sin_p::Float64,E0::Float64,g0::Int64,Σt::Matrix{Float64},Tmt::Matrix{Float64},Sb::Matrix{Float64},Eb::Vector{Float64},ΔE::Vector{Float64},widths::Vector{Vector{Float64}},mat3::Array{Int64,3},Np::Int64,Ndims::Int64,r0::Vector{Float64},Ω::Vector{Float64},bounds::Vector{Vector{Float64}})
    Ng = length(ΔE); Lmax = maximum(pl)
    fac = zeros(Lmax+1)
    g = g0; E = E0; T_att = 1.0; τ_T = 0.0
    z̄ = 0.0; Z2 = 0.0; X2 = 0.0; Jz = 0.0; Jx = 0.0; g_pos = 0.0
    iseg = 1
    if isempty(segs) return end
    ix,iy,iz,ℓseg = segs[iseg]; d_geom = ℓseg
    while g ≤ Ng && T_att > 1e-300
        m = mat3[ix,iy,iz]
        E_lo = Eb[g+1]
        k = (Sb[g,m]-Sb[g+1,m])/(Eb[g]-E_lo)
        S_E = Sb[g+1,m]+k*(E-E_lo)
        Δs_grp = abs(k)*ΔE[g] > 1e-12*S_E ? log(S_E/Sb[g+1,m])/k : (E-E_lo)/S_E
        T_gm = Tmt[g,m]
        if T_gm > 0
            # Invert the D₁-compressed geometric advance to reach the segment end:
            # d_geom = sin_p·e^{-2τ_T}·(1-e^{-2TΔs})/(2T).
            arg = 1.0-2*T_gm*(d_geom/sin_p)*exp(2*τ_T)
            Δs_seg = arg > 0 ? -log(arg)/(2*T_gm) : Inf
        else
            Δs_seg = d_geom/sin_p
        end
        crosses_group = Δs_grp ≤ Δs_seg
        Δs = crosses_group ? Δs_grp : Δs_seg
        Σ̄ = Σt[g,m]
        for l in range(0,Lmax)
            Σ_eff = Σ̄+l*(l+1)*T_gm
            fac[l+1] = exp(-l*(l+1)*τ_T)*(Σ_eff*Δs > 1e-12 ? (1-exp(-Σ_eff*Δs))/Σ_eff : Δs)
        end
        Δd = T_gm > 0 ? sin_p*exp(-2*τ_T)*(1-exp(-2*T_gm*Δs))/(2*T_gm) : sin_p*exp(-2*τ_T)*Δs
        τ_T,z̄,Z2,X2,Jz,Jx = fcs_lewis_update(τ_T,z̄,Z2,X2,Jz,Jx,T_gm,Δs)
        # Deposit at the substep geometric midpoint, spread by the Lewis kernel.
        jlo,wsp = fcs_splat_windows(g_pos+Δd/2,max(Z2-z̄^2,0.0),max(X2,0.0),r0,Ω,sin_p,bounds,(ix,iy,iz),Ndims)
        for (kz,wkz) in enumerate(wsp[3]), (ky,wky) in enumerate(wsp[2]), (kx,wkx) in enumerate(wsp[1])
            K = wkx*wky*wkz
            if K ≤ 1e-300 continue end
            jx = jlo[1]+kx-1; jy = jlo[2]+ky-1; jz = jlo[3]+kz-1
            Vc = widths[1][jx]*(Ndims ≥ 2 ? widths[2][jy] : 1.0)*(Ndims ≥ 3 ? widths[3][jz] : 1.0)
            w_seg = (wA*K/Vc)*T_att
            for p in range(1,Np)
                if 𝒴[p] != 0.0 φ_u[g,p,jx,jy,jz] += w_seg*𝒴[p]*fac[pl[p]+1] end
            end
        end
        T_att *= exp(-Σ̄*Δs)
        g_pos += Δd
        if crosses_group
            d_geom -= Δd
            g += 1; E = E_lo
            if g > Ng
                jlo,wsp = fcs_splat_windows(g_pos,max(Z2-z̄^2,0.0),max(X2,0.0),r0,Ω,sin_p,bounds,(ix,iy,iz),Ndims)
                for (kz,wkz) in enumerate(wsp[3]), (ky,wky) in enumerate(wsp[2]), (kx,wkx) in enumerate(wsp[1])
                    K = wkx*wky*wkz
                    if K ≤ 1e-300 continue end
                    jx = jlo[1]+kx-1; jy = jlo[2]+ky-1; jz = jlo[3]+kz-1
                    Vc = widths[1][jx]*(Ndims ≥ 2 ? widths[2][jy] : 1.0)*(Ndims ≥ 3 ? widths[3][jz] : 1.0)
                    w_cut = (wA*K/Vc)*T_att/Sb[Ng+1,m]*ΔE[Ng]
                    for p in range(1,Np)
                        if 𝒴[p] != 0.0 φ_cut[p,jx,jy,jz] += w_cut*𝒴[p]*exp(-pl[p]*(pl[p]+1)*τ_T) end
                    end
                end
                return
            end
        else
            if abs(k)*ΔE[g] > 1e-12*S_E
                E = E_lo+(S_E*exp(-k*Δs)-Sb[g+1,m])/k
            else
                E -= S_E*Δs
            end
            iseg += 1
            if iseg > length(segs) return end
            ix,iy,iz,ℓseg = segs[iseg]; d_geom = ℓseg
        end
    end
end

"""
    fcs_lewis_update(τ_T,z̄,Z2,X2,Jz,Jx,T,Δs)

Advance the exact Lewis moments of a Goudsmit--Saunderson pencil over a substep of
true pathlength `Δs` with constant restricted momentum transfer `T`: mean forward
displacement `z̄ = ∫D₁ds` along the beam axis, second moments `Z2 = ⟨z²⟩` (beam
axis) and `X2 = ⟨x²⟩` (per transverse axis) via their correlation integrals
`Jz,Jx` (direction correlations propagate exactly by `D̃₁ = e^(-2∫Tds)`, with
`⟨μ²⟩ = (1+2D₂)/3`, `⟨u_x²⟩ = (1-D₂)/3`, `D₂ = e^(-6τ_T)`). Closed form per
substep; frozen-coefficient limit for `TΔs ≪ 1`.

# Input Argument(s)
- `τ_T,z̄,Z2,X2,Jz,Jx::Float64` : Lewis state at the substep entry.
- `T::Float64` : restricted momentum transfer (Radiant convention `λ_ℓ = ℓ(ℓ+1)T`).
- `Δs::Float64` : true pathlength of the substep.

# Output Argument(s)
- `τ_T,z̄,Z2,X2,Jz,Jx::Float64` : Lewis state at the substep exit.

"""
function fcs_lewis_update(τ_T::Float64,z̄::Float64,Z2::Float64,X2::Float64,Jz::Float64,Jx::Float64,T::Float64,Δs::Float64)
    A = exp(-2*τ_T); B = exp(-6*τ_T)
    if T*Δs > 1e-8
        q2 = exp(-2*T*Δs); q6 = exp(-6*T*Δs)
        E2 = (1-q2)/(2*T); E6 = (1-q6)/(6*T); Q4 = (q2-q6)/(4*T)
        z̄ += A*E2
        Z2 += 2*Jz*E2 + (Δs-E2)/(3*T) + B/(3*T)*(E2-E6)
        X2 += 2*Jx*E2 + (Δs-E2)/(3*T) - B/(6*T)*(E2-E6)
        Jz = Jz*q2 + (E2+2*B*Q4)/3
        Jx = Jx*q2 + (E2-B*Q4)/3
    else
        z̄ += A*Δs
        Z2 += 2*Jz*Δs + (1+2*B)/3*Δs^2
        X2 += 2*Jx*Δs + (1-B)/3*Δs^2
        Jz += (1+2*B)/3*Δs
        Jx += (1-B)/3*Δs
    end
    return τ_T+T*Δs,z̄,Z2,X2,Jz,Jx
end

"""
    fcs_splat_windows(g_mid,σz2,σx2,r0,Ω,sin_p,bounds,host,Ndims)

Per-axis Gaussian kernel windows for a Lewis-spread deposit centered at geometric
ray distance `g_mid`: the kernel covariance is projected on each grid axis
(`σ_a² = σ⊥²+(σ_z²-σ⊥²)·Ω_a²`, the exact marginal variances; the cross-covariance
of oblique rays is neglected) and integrated per voxel by `fcs_splat_axis`.

# Input Argument(s)
- `g_mid::Float64` : geometric distance of the kernel center along the ray.
- `σz2,σx2::Float64` : longitudinal and transverse Lewis variances.
- `r0::Vector{Float64}`, `Ω::Vector{Float64}`, `sin_p::Float64` : ray geometry.
- `bounds::Vector{Vector{Float64}}` : voxel boundaries per active axis.
- `host::Tuple` : walker voxel (inactive axes and degenerate kernels).
- `Ndims::Int64` : number of active spatial dimensions.

# Output Argument(s)
- `jlo::Vector{Int64}` : first voxel of the window per axis.
- `w::Vector{Vector{Float64}}` : kernel mass per voxel of the window per axis.

"""
function fcs_splat_windows(g_mid::Float64,σz2::Float64,σx2::Float64,r0::Vector{Float64},Ω::Vector{Float64},sin_p::Float64,bounds::Vector{Vector{Float64}},host::Tuple{Int64,Int64,Int64},Ndims::Int64)
    jlo = Vector{Int64}(undef,3); w = Vector{Vector{Float64}}(undef,3)
    for a in range(1,3)
        if a ≤ Ndims
            σ_a = sqrt(max(σx2+(σz2-σx2)*Ω[a]^2,0.0))
            jlo[a],w[a] = fcs_splat_axis(r0[a]+Ω[a]/sin_p*g_mid,σ_a,bounds[a],host[a])
        else
            jlo[a] = host[a]; w[a] = [1.0]
        end
    end
    return jlo,w
end

"""
    fcs_splat_axis(r̄::Float64,σ::Float64,b::Vector{Float64},jhost::Int64)

Gaussian kernel weights along one grid axis: mass of `N(r̄,σ²)` integrated over
each voxel of `b` within `r̄±6σ` (mass beyond the mesh is leaked, consistent with
void boundaries). A kernel much thinner than the host voxel collapses to it.

# Input Argument(s)
- `r̄::Float64`, `σ::Float64` : kernel center and standard deviation.
- `b::Vector{Float64}` : voxel boundaries along the axis.
- `jhost::Int64` : voxel containing the deposit (degenerate-kernel fallback).

# Output Argument(s)
- `jlo::Int64` : first voxel of the window.
- `w::Vector{Float64}` : kernel mass per voxel of the window.

"""
function fcs_splat_axis(r̄::Float64,σ::Float64,b::Vector{Float64},jhost::Int64)
    if 6*σ ≤ 0.01*(b[jhost+1]-b[jhost]) return jhost,[1.0] end
    N = length(b)-1
    jlo = clamp(searchsortedlast(b,r̄-6*σ),1,N)
    jhi = clamp(searchsortedlast(b,r̄+6*σ),1,N)
    w = Vector{Float64}(undef,jhi-jlo+1)
    c = 1/(sqrt(2)*σ)
    Φ = 0.5*erf((b[jlo]-r̄)*c)
    for j in range(jlo,jhi)
        Φn = 0.5*erf((b[j+1]-r̄)*c)
        w[j-jlo+1] = Φn-Φ
        Φ = Φn
    end
    return jlo,w
end
