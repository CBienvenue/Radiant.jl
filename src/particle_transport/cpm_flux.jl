"""
    compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::CPM,
    source::Source,electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field())

Solve the transport equation using the collision-probability method (CPM) for a given
particle, in one-dimensional Cartesian geometry.

The azimuthally-symmetric (`m = 0`) volume flux is expanded in Legendre polynomials up to
degree `Lp = solver.get_legendre_order()`. For each energy group the volume--volume
collision-probability matrix ``P_{vv}`` is assembled (see
[`cpm_collision_matrix`](@ref)); the in-group scattering closure
``(\\mathbb{I} - P_{vv}\\Sigma_s)\\vec{\\phi} = P_{vv}\\vec{Q}`` is solved directly, where the
source ``\\vec{Q}`` gathers the fixed volume source and the out-of-group scattering source.
Vacuum boundary conditions are assumed (the incoming partial currents vanish, so the surface
probability matrices do not contribute).

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross-section informations.
- `geometry::Geometry` : geometry informations.
- `solver::CPM` : collision-probability method informations.
- `source::Source` : source informations.
- `electromagnetic_field::Electromagnetic_Field` : external electromagnetic field (optional,
  unused by the CPM solver).

# Output Argument(s)
- `flux::Flux_Per_Particle` : flux informations.

# Reference(s)
- TR-02, *Collision Probability Methods* (Radiant technical report).

"""
function compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::CPM,source::Source,electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field())

#----
# Geometry data
#----
Ndims = geometry.get_dimension()
geo_type = geometry.get_type()
if geo_type != "cartesian" error("The CPM solver is only available in Cartesian geometry.") end
if Ndims != 1 error("The CPM solver is currently only available in 1D Cartesian geometry.") end
Ns = geometry.get_number_of_voxels()
Δs = geometry.get_voxels_width()
mat = geometry.get_material_per_voxel()
boundary_conditions = geometry.get_boundary_conditions()
if any(boundary_conditions .== 2) error("The CPM solver does not support periodic boundary conditions.") end
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
    # P̃vv = Pvv + Pvs(I-A Pss)⁻¹A Psv and the LU factorization of the in-group scattering operator
    # (I - P̃vv Σs_in). Both are independent of the outer iteration, so this is done once.
    Ptvv = Vector{Matrix{Float64}}(undef,Ng)
    Fact = Vector{Any}(undef,Ng)
    for ig in range(1,Ng)
        Σcell = [Σtot[ig,mat[i,1,1]] for i in range(1,Nx)]
        Pvv = cpm_collision_matrix(Σcell,Δx,Lp)
        if is_reflective
            Pvs,Psv,Pss = cpm_surface_matrices(Σcell,Δx,Lp,Nν)
            Sv = Nν+1
            A = Diagonal(vcat(fill(βL,Sv),fill(βR,Sv)))
            Pvv = Pvv .+ Pvs*((Matrix{Float64}(I,2Sv,2Sv) .- A*Pss) \ (A*Psv))
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
            φ = Fact[ig] \ (Ptvv[ig] * Q)
            for i in range(1,Nx), p in range(1,P)
                𝚽l[ig,p,1,i,1,1] = φ[idx(i,p)]
            end
        end
        has_upscatter || break
        if norm(𝚽l .- 𝚽l⁻) / max(norm(𝚽l),eps()) < ϵ_out break end
    end

else  # mode == "sweeping"

    Sv = Nν+1
    # Per (group, cell) single-voxel response blocks (source/incoming ↦ flux/outgoing).
    cells = [[cpm_cell_response(Σtot[ig,mat[i,1,1]],Δx[i],Lp,Nν) for i in range(1,Nx)] for ig in range(1,Ng)]

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
            φ,ρ_in[ig] = cpm_sweep_1D(Qext,cells[ig],Σs_in,βL,βR,P,Sv,Nx,𝒜,I_max,ϵ_out,ig,gmres_restart,anderson_depth)
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
