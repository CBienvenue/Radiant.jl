"""
    compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::CP,
    source::Source,electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field())

Solve the transport equation using the collision-probability method (CP) for a given
particle, in one-dimensional Cartesian geometry.

The azimuthally-symmetric (`m = 0`) volume flux is expanded in Legendre polynomials up to
degree `Lp = solver.get_legendre_order()`. For each energy group the volume--volume
collision-probability matrix ``P_{vv}`` is assembled (see
[`cp_collision_matrix`](@ref)); the in-group scattering closure
``(\\mathbb{I} - P_{vv}\\Sigma_s)\\vec{\\phi} = P_{vv}\\vec{Q}`` is solved directly, where the
source ``\\vec{Q}`` gathers the fixed volume source and the out-of-group scattering source.

Vacuum, reflective and surface-source boundaries are supported. A fixed incoming boundary flux
(surface source) enters through its half-range moments ``\\vec{g}`` in the surface basis
``\\bar{R}_\\ell``: it adds the term ``b = P_{vs}(\\mathbb{I} - A\\,P_{ss})^{-1}\\vec{g}`` to the
right-hand side (``A`` the albedo of the reflective faces), or, in the sweeping mode, seeds the
incoming interface currents. With vacuum boundaries and no surface source the incoming partial
currents vanish and the surface probability matrices do not contribute.

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross-section informations.
- `geometry::Geometry` : geometry informations.
- `solver::CP` : collision-probability method informations.
- `source::Source` : source informations.
- `electromagnetic_field::Electromagnetic_Field` : external electromagnetic field (optional,
  unused by the CP solver).

# Output Argument(s)
- `flux::Flux_Per_Particle` : flux informations.

# Reference(s)
- TR-02, *Collision Probability Methods* (Radiant technical report).

"""
function compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::CP,source::Source,electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field())

#----
# Geometry data
#----
Ndims = geometry.get_dimension()
geo_type = geometry.get_type()
if geo_type != "cartesian" error("The CP solver is only available in Cartesian geometry.") end
if Ndims != 1 error("The CP solver is currently only available in 1D Cartesian geometry.") end
Ns = geometry.get_number_of_voxels()
Δs = geometry.get_voxels_width()
mat = geometry.get_material_per_voxel()
boundary_conditions = geometry.get_boundary_conditions()
if any(boundary_conditions .== 2) error("The CP solver does not support periodic boundary conditions.") end
Nx = Ns[1]
Δx = Δs[1]

# Albedo (reflection) coefficients per face : 0 = vacuum, 1 = reflective.
βL = boundary_conditions[1] == 1 ? 1.0 : 0.0   # x- (left) face
βR = boundary_conditions[2] == 1 ? 1.0 : 0.0   # x+ (right) face
is_reflective = (βL != 0.0) || (βR != 0.0)

#----
# Solver / angular basis
#----
part = solver.get_particle()
solver.get_solver_type()   # validates that the solver type is supported (BTE)
Lp = solver.get_legendre_order()
P = Lp+1
Nν = solver.get_surface_order()

#----
# Cross sections
#----
Ng = cross_sections.get_number_of_groups(part)
Nmat = cross_sections.get_number_of_materials()
Σtot = cross_sections.get_total(part)                 # [Ng,Nmat]
Σs = cross_sections.get_scattering(part,part,Lp)      # [Nmat,Ng,Ng,Lp+1] (from group, to group, ℓ)

println(">>>Particle: $(get_type(part)) <<<")

#----
# Fixed sources
#----
volume_sources = source.get_volume_sources()          # [Ng,P,1,Nx,1,1]
surface_sources = source.get_surface_sources()        # [Ng,Nνsrc+1,2] : half-range moments per face

# Incoming boundary flux moments in the half-range basis R̄_ℓ, stacked as
# [left face ℓ = 0…Nν, right face ℓ = 0…Nν] (S = Nν+1 moments per face). The surface source
# supplies its moments for ℓ ≤ min(Nν, Nνsrc); the remaining ones stay zero.
Sv = Nν+1
Lsrc = size(surface_sources,2)-1
has_surface_source = false
for ig in range(1,Ng), l in range(1,Lsrc+1), f in range(1,2)
    if surface_sources[ig,l,f] != 0.0 has_surface_source = true end
end
boundary_vector(ig) = begin
    g = zeros(2*Sv)
    for l in range(0,min(Nν,Lsrc))
        g[l+1]    = convert(Float64,surface_sources[ig,l+1,1])   # left face  (X-)
        g[Sv+l+1] = convert(Float64,surface_sources[ig,l+1,2])   # right face (X+)
    end
    # The surface-source moments are the angular-flux half-range moments J_ℓ = ∫ψ R̄_ℓ dμ̂
    # (Radiant's convention, shared with SN). The CP boundary-coupling matrices Pvs/Css consume
    # the azimuthally-integrated interface unknown ∫∫ψ R̄_ℓ dμ̂ dφ = 2π J_ℓ, so divide by 2π.
    return g ./ (2π)
end

#----
# Flux calculations
#----
𝚽l = zeros(Ng,P,1,Ns[1],Ns[2],Ns[3])
ρ_in = -ones(Ng)
mode = solver.get_mode()
𝒜 = solver.get_acceleration()
I_max = solver.get_maximum_iteration()
ϵ_out = solver.get_convergence_criterion()
gmres_restart = solver.get_gmres_restart()
anderson_depth = solver.get_anderson_depth()

# Index of moment p (1-based) in region i (1-based) within the flat (Nx*P) vectors.
idx(i,p) = (i-1)*P + p

# Up-scattering detection : a nonzero transfer from a higher-index group into a lower-index one.
# Down-scattering is resolved in a single ordered group pass; up-scattering additionally requires
# an outer Gauss-Seidel iteration over the energy groups.
has_upscatter = false
for ig in range(1,Ng), gi in range(ig+1,Ng), n in range(1,Nmat)
    if Σs[n,gi,ig,1] != 0.0 has_upscatter = true end
end
Nout = has_upscatter ? I_max : 1

if mode == "global"

    # Per energy group, precompute the boundary-closed collision matrix
    # P̃vv = Pvv + Pvs(I-A Pss)⁻¹A Psv, the LU factorization of the in-group scattering operator
    # (I - P̃vv Σs_in), and the fixed surface-source contribution to the volume flux
    # b = Pvs(I-A Pss)⁻¹g. All three are independent of the outer iteration, so this is done once.
    Ptvv = Vector{Matrix{Float64}}(undef,Ng)
    Fact = Vector{Any}(undef,Ng)
    bsrc = [zeros(Nx*P) for _ in range(1,Ng)]
    for ig in range(1,Ng)
        Σcell = [Σtot[ig,mat[i,1,1]] for i in range(1,Nx)]
        Pvv = cp_collision_matrix(Σcell,Δx,Lp)
        if is_reflective || has_surface_source
            Pvs,Psv,Pss = cp_surface_matrices(Σcell,Δx,Lp,Nν)
            A = Diagonal(vcat(fill(βL,Sv),fill(βR,Sv)))
            M = Matrix{Float64}(I,2Sv,2Sv) .- A*Pss          # (I − A Pss)
            if is_reflective      Pvv = Pvv .+ Pvs*(M \ (A*Psv)) end
            if has_surface_source bsrc[ig] = Pvs*(M \ boundary_vector(ig)) end
        end
        Σs_in = [Σs[mat[j,1,1],ig,ig,q] for j in range(1,Nx) for q in range(1,P)]
        Ptvv[ig] = Pvv
        Fact[ig] = lu(Matrix{Float64}(I,Nx*P,Nx*P) .- Pvv * Diagonal(Σs_in))
    end

    for i_out in range(1,Nout)
        𝚽l⁻ = copy(𝚽l)
        for ig in range(1,Ng)
            # Source moments : fixed volume source + out-of-group scattering (current fluxes).
            Q = zeros(Nx*P)
            for j in range(1,Nx), q in range(1,P)
                m = mat[j,1,1]
                Q[idx(j,q)] += volume_sources[ig,q,1,j,1,1]
                for gi in range(1,Ng)
                    if gi != ig
                        Q[idx(j,q)] += Σs[m,gi,ig,q] * 𝚽l[gi,q,1,j,1,1]
                    end
                end
            end
            φ = Fact[ig] \ (Ptvv[ig] * Q .+ bsrc[ig])
            for i in range(1,Nx), p in range(1,P)
                𝚽l[ig,p,1,i,1,1] = φ[idx(i,p)]
            end
        end
        has_upscatter || break
        if norm(𝚽l .- 𝚽l⁻) / max(norm(𝚽l),eps()) < ϵ_out break end
    end

else  # mode == "sweeping"

    # Per (group, cell) single-voxel response blocks (source/incoming ↦ flux/outgoing).
    cells = [[cp_cell_response(Σtot[ig,mat[i,1,1]],Δx[i],Lp,Nν) for i in range(1,Nx)] for ig in range(1,Ng)]

    for i_out in range(1,Nout)
        𝚽l⁻ = copy(𝚽l)
        for ig in range(1,Ng)
            # External source moments (P × Nx) : fixed volume source + out-of-group scattering.
            Qext = zeros(P,Nx)
            for j in range(1,Nx), q in range(1,P)
                m = mat[j,1,1]
                Qext[q,j] += volume_sources[ig,q,1,j,1,1]
                for gi in range(1,Ng)
                    if gi != ig
                        Qext[q,j] += Σs[m,gi,ig,q] * 𝚽l[gi,q,1,j,1,1]
                    end
                end
            end
            # In-group scattering Legendre moments (P × Nx).
            Σs_in = [Σs[mat[j,1,1],ig,ig,q] for q in range(1,P), j in range(1,Nx)]
            # Fixed incoming boundary current (left / right) from the surface source.
            g = boundary_vector(ig)
            φ,ρ_in[ig] = cp_sweep_1D(Qext,cells[ig],Σs_in,βL,βR,g[1:Sv],g[Sv+1:2Sv],P,Sv,Nx,𝒜,I_max,ϵ_out,ig,gmres_restart,anderson_depth)
            for i in range(1,Nx), p in range(1,P)
                𝚽l[ig,p,1,i,1,1] = φ[p,i]
            end
        end
        has_upscatter || break
        if norm(𝚽l .- 𝚽l⁻) / max(norm(𝚽l),eps()) < ϵ_out break end
    end

end

#----
# Save flux
#----
flux = Flux_Per_Particle(part)
flux.add_flux(𝚽l)
flux.add_spectral_radius(ρ_in)

return flux

end
