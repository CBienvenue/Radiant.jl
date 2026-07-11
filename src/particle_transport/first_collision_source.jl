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
- `ss::Surface_Source` : surface source with `beam_treatment = "first-collision"`.

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
    if lowercase(geometry.get_type()) != "cartesian" || geometry.get_dimension() != 1
        error("The first-collision source treatment is only available in 1D cartesian geometries.")
    end
    bc = geometry.get_boundary_conditions()
    if bc[1] != 0 || bc[2] != 0
        error("The first-collision source treatment requires void x-boundaries (the analytic uncollided flux is neither reflected nor wrapped).")
    end
    solver_type,is_CSD = solver.get_solver_type()
    if solver_type ∉ [1,2,3]
        error("The first-collision source treatment is only available for the BTE, BFP and BCSD solvers.")
    end
    if ismissing(ss.energy_group) error("The energy group of the surface source is not defined.") end
    if ismissing(ss.location) error("The location of the surface source is not defined.") end
    face = uppercase(ss.location)
    if face ∉ ["X-","X+"]
        error("The first-collision source treatment is only available on the X- and X+ faces.")
    end
    is_delta = ~ismissing(ss.direction)
    if ~is_delta && ismissing(ss.angular_moments)
        error("A first-collision surface source requires either a direction (set_direction) or half-range angular moments (set_angular_moments).")
    end
    if is_delta
        μs = ss.direction[1]
        if (face == "X-" && μs < 1e-12) || (face == "X+" && μs > -1e-12)
            error("The monodirectional surface source must be incoming (and not grazing) with respect to its face.")
        end
    end

    #----
    # Angular basis, geometry and material data
    #----
    Np,pl,pm,is_SPH = fcs_angular_basis(solver,1)
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
            𝒴 = fcs_direction_basis(is_SPH,pl,pm,ss.direction)
            fcs_bte_delta!(φ_u,𝒴,I,abs(ss.direction[1]),face,g0,Σt[g0,:],Δx,mat,Nx,Np)
        else
            fcs_bte_moments!(φ_u,ss.angular_moments,I,face,g0,pl,pm,is_SPH,Σt[g0,:],Δx,mat,Nx,Np)
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
            𝒴 = fcs_direction_basis(is_SPH,pl,pm,ss.direction)
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
                𝒴 = fcs_column_basis(is_SPH,pl,pm,face == "X-" ? μ̂ : -μ̂)
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
        if Qdims != 1 error("The first-collision source with the SN solver requires a 1D quadrature.") end
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
    fcs_direction_basis(is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},Ω::Vector{Float64})

Full-range moment values `𝒴_p` of the angular δ-function at direction `Ω`:
`P_ℓ(μₛ)` (Legendre/SN bases) or `R_ℓm(μₛ,ϕₛ)` (spherical harmonics).

# Input Argument(s)
- `is_SPH::Bool` : spherical-harmonics (true) or Legendre-type (false) basis.
- `pl::Vector{Int64}`, `pm::Vector{Int64}` : moment orders.
- `Ω::Vector{Float64}` : direction cosines [μ,η,ξ].

# Output Argument(s)
- `𝒴::Vector{Float64}` : basis values per moment.

"""
function fcs_direction_basis(is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},Ω::Vector{Float64})
    Np = length(pl)
    Lmax = maximum(pl)
    𝒴 = zeros(Np)
    if is_SPH
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
    fcs_column_basis(is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},μ::Float64)

Full-range moment values `P_ℓ(μ)` of an azimuthally-symmetric column at direction
cosine `μ`; only the m = 0 slots are filled in the spherical-harmonics basis.

# Input Argument(s)
- `is_SPH::Bool` : spherical-harmonics (true) or Legendre-type (false) basis.
- `pl::Vector{Int64}`, `pm::Vector{Int64}` : moment orders.
- `μ::Float64` : direction cosine.

# Output Argument(s)
- `𝒴::Vector{Float64}` : basis values per moment.

"""
function fcs_column_basis(is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},μ::Float64)
    Np = length(pl)
    𝒴 = zeros(Np)
    Pls = legendre_polynomials_up_to_L(maximum(pl),μ)
    for p in range(1,Np)
        if ~is_SPH || pm[p] == 0
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
    fcs_bte_moments!(φ_u,g_mom,I,face,g0,pl,pm,is_SPH,Σt_g,Δx,mat,Nx,Np)

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
- `pl::Vector{Int64}`, `pm::Vector{Int64}`, `is_SPH::Bool` : moment basis.
- `Σt_g::Vector{Float64}` : total cross section per material in the emission group.
- `Δx::Vector{Float64}`, `mat::Vector{Int64}`, `Nx::Int64`, `Np::Int64` : mesh data.

# Output Argument(s)
N/A

"""
function fcs_bte_moments!(φ_u::Array{Float64,3},g_mom::Vector{Float64},I::Float64,face::String,g0::Int64,pl::Vector{Int64},pm::Vector{Int64},is_SPH::Bool,Σt_g::Vector{Float64},Δx::Vector{Float64},mat::Vector{Int64},Nx::Int64,Np::Int64)
    Lmax = maximum(pl)
    L_src = length(g_mom)-1
    # Moment slots receiving the order-ℓ value (m = 0 only in the SH basis).
    slots = [Int64[] for l in 0:Lmax]
    for p in range(1,Np)
        if ~is_SPH || pm[p] == 0 push!(slots[pl[p]+1],p) end
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
