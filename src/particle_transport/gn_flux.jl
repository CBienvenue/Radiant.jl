"""
    compute_flux(cross_sections::Cross_Sections,geometry::Geometry,
    solver::GN,source::Source)

Solve the transport equation using the spherical Galerkin harmonics (GN) method for a given particle.  

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross section informations.
- `geometry::Geometry` : geometry informations.
- `solver::GN` : spherical Galerkin harmonics informations.
- `source::Source` : source informations.
- `electromagnetic_field::Electromagnetic_Field` : external electromagnetic field (optional,
  defaults to no field; external fields are not supported with the GN solver).

# Output Argument(s)
- `flux::Flux_Per_Particle`: flux informations.

# Reference(s)
N/A

"""
function compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::GN,source::Source,electromagnetic_field::Electromagnetic_Field=Electromagnetic_Field())

# External electromagnetic fields are only implemented for the SN solver.
if electromagnetic_field.get_electromagnetic_field()[1] error("External electromagnetic fields are not supported with the GN solver.") end

#----
# Geometry data
#----
Ndims = geometry.get_dimension()
geo_type = geometry.get_type()
if geo_type != "cartesian" error("Transport of particles in",geo_type," is unavailable.") end
Ns = geometry.get_number_of_voxels()
Δs = geometry.get_voxels_width()
mat = geometry.get_material_per_voxel()
boundary_conditions = geometry.get_boundary_conditions()

#----
# Preparation of angular discretisation
#----
L = solver.get_legendre_order()
L_elem = solver.get_legendre_order_local()
Nv = solver.get_subdivision()
tiling = solver.get_tiling()
# `z_fold` folds the angular domain by the symmetry of the reduced-dimension problem:
# in 1D (spherical harmonics) the azimuthally-symmetric flux collapses onto two
# half-spheres; in 2D the μ_z → -μ_z symmetry collapses the sphere onto four
# quadrants (with Nv = 1). With `z_fold = false` the full octant tiling is used in
# every dimension (the general GN treatment, with azimuthal subdivision).
fold    = solver.get_z_fold()
polynomial_basis = solver.get_polynomial_basis(Ndims)
fold_1D = (Ndims == 1 && fold && polynomial_basis == "spherical-harmonics")  # 1D SH → two half-spheres
z_fold  = (Ndims == 2 && Nv == 1 && fold)            # 2D → four quadrants
if polynomial_basis == "legendre"
    if Ndims != 1 error("Legendre basis is only available in 1D.") end
    is_SPH = false
    Np,Nq,Mll = patch_to_full_range_matrix_legendre(L,L_elem,Nv)
    pl = collect(0:L)
    pm = zeros(Int64,Np)
    𝒩 = gn_weights_legendre_1D(L_elem,Nv)
elseif polynomial_basis == "spherical-harmonics"
    is_SPH = true
    if fold_1D
        # Two half-spheres (azimuthally collapsed), Nv polar bands per hemisphere.
        Np,Nq,Mll = patch_to_full_range_matrix_spherical_harmonics_1D(L,L_elem,Nv)
        pl,pm = spherical_harmonics_indices(L)
        𝒩 = gn_weights_spherical_harmonics_1D(L_elem,Nv)
    elseif z_fold
        # Four quadrants (z-symmetry fold).
        Np,Nq,Mll = patch_to_full_range_matrix_spherical_harmonics_2D_quarter(L,L_elem)
        pl,pm = spherical_harmonics_indices(L)
        𝒩 = gn_weights_spherical_harmonics_2D_quarter(L_elem)
    else
        # Full octant tiling with azimuthal subdivision (general GN, any dimension).
        Np,Nq,Mll = patch_to_full_range_matrix_spherical_harmonics(L,L_elem,Nv;tiling=tiling)
        pl,pm = spherical_harmonics_indices(L)
        𝒩 = gn_weights_spherical_harmonics(L_elem,Nv,Ndims;tiling=tiling)
    end
else
    error("Unknown polynomial basis.")
end

#----
# Preparation of cross sections
#----

part = solver.get_particle()
solver_type,is_CSD = solver.get_solver_type()
Nmat = cross_sections.get_number_of_materials()
Ng = cross_sections.get_number_of_groups(part)
if is_CSD
    ΔE = cross_sections.get_energy_width(part)
    E = cross_sections.get_energies(part)
    Eb = cross_sections.get_energy_boundaries(part)
end

isFC = solver.get_is_full_coupling()
schemes,𝒪,Nm = solver.get_schemes(geometry,isFC)
ω,𝒞,𝒲 = scheme_weights_gn(𝒪,schemes,Ndims,is_CSD)

println(">>>Particle: $(get_type(part)) <<<")

# Total cross sections
Σtot = zeros(Ng,Nmat)
if solver_type ∈ [4,5] 
    Σtot = cross_sections.get_absorption(part)
else
    Σtot = cross_sections.get_total(part)
end

# Scattering cross sections
Σs = zeros(Nmat,Ng,Ng,L+1)
if solver_type ∉ [4,5]
    Σs = cross_sections.get_scattering(part,part,L)
end

# Stopping powers
if is_CSD
    S⁻ = zeros(Ng,Nmat); S⁺ = zeros(Ng,Nmat)
    Sb = cross_sections.get_boundary_stopping_powers(part)
    for n in range(1,Nmat)
        S⁻[:,n] = Sb[1:Ng,n] ; S⁺[:,n] = Sb[2:Ng+1,n]
    end
    S = zeros(Ng,Nmat,𝒪[4])
    for n in range(1,Nmat), ig in range(1,Ng)
        S[ig,n,1] = (S⁻[ig,n]+S⁺[ig,n])/2
        if (𝒪[4] > 1) S[ig,n,2] = (S⁻[ig,n]-S⁺[ig,n])/(2*sqrt(3)) end
    end
end

# Momentum transfer
if solver_type ∈ [2,4]
    T = zeros(Ng,Nmat)
    T = cross_sections.get_momentum_transfer(part)
    fokker_planck_type = solver.get_angular_fokker_planck()
    if fold_1D
        # 1D spherical-harmonics over two half-spheres: only the dimension-agnostic
        # galerkin operator is supported (the finite-difference patch graph assumes
        # the octant tiling). This is the same capability as the DPN solver.
        if fokker_planck_type != "galerkin" error("The folded 1D spherical-harmonics GN basis (two half-spheres) only supports the \"galerkin\" angular Fokker-Planck operator; use set_z_fold(false) for the octant tiling.") end
        ℳ,λ₀ = fokker_planck_scattering_matrix("galerkin",pl,Np)
    elseif is_SPH
        ℳ,λ₀ = fokker_planck_scattering_matrix(fokker_planck_type,pl,Np;L=L,L_elem=L_elem,Nv=Nv,Ndims=Ndims,tiling=tiling,Mll=Mll)
    elseif fokker_planck_type == "finite-difference"
        ℳ,λ₀ = fokker_planck_finite_difference_gn_legendre_1D(L,L_elem,Nv,Mll)
    else
        ℳ,λ₀ = fokker_planck_scattering_matrix("galerkin",pl,Np)
    end
    Σtot .+= T .* λ₀
end

# Elastic-free approximation
if solver_type == 6
    for n in range(1,Nmat), ig in range(1,Ng)
        Σtot[ig,n] -= Σs[n,ig,ig,1]
    end
end

#----
# Acceleration solver
#----

𝒜 = solver.get_acceleration()
gmres_restart = solver.get_gmres_restart()
anderson_depth = solver.get_anderson_depth()

#----
# Fixed sources
#----

L_surf = 15
if is_SPH
    Np_surf = spherical_harmonics_number_basis(L_surf)
else
    Np_surf = L_surf + 1
end
surface_sources = source.get_surface_sources()
volume_sources = source.get_volume_sources()
Np_source = Int64(min(Np_surf,length(surface_sources[1,:,1])))
if fold_1D
    Mll_surf = patch_to_half_range_matrix_spherical_harmonics_1D(L_surf,L_elem,Nv)
elseif z_fold
    Mll_surf = patch_to_half_range_matrix_spherical_harmonics_2D_quarter(L_surf,L_elem)
elseif is_SPH
    Mll_surf = patch_to_half_range_matrix_spherical_harmonics(L_surf,L_elem,Nv,Ndims;tiling=tiling)
else
    Mll_surf = patch_to_half_range_matrix_legendre(L_surf,L_elem,Nv,Ndims)
end

#----
# Boundary conditions
#----
if any(x->x == 1,boundary_conditions) # If reflective, construct the reflection matrices
    if is_SPH
        Rpq = pos_to_neg_half_range_matrix_spherical_harmonics(L_surf,Ndims)
    else
        Rpq = pos_to_neg_half_range_matrix_legendre(L_surf,Ndims)
    end
else
    Rpq = Array{Float64}(undef,0,0,0)
end

#----
# Flux calculations
#----

ϵ_max = solver.get_convergence_criterion()
I_max = solver.get_maximum_iteration()

# Initialization flux
𝚽l = zeros(Ng,Np,Nm[5],Ns[1],Ns[2],Ns[3])
if is_CSD 𝚽cutoff = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3]) end

# All-group iteration
i_out = 1
is_outer_convergence = false
ϵ_out = Inf
is_outer_iteration = false
Ntot = 0
ρ_in = -ones(Ng) # In-group spectral radius (per energy group, last outer iteration)
if is_outer_iteration 𝚽l⁻ = zeros(Ng,Np,Nm[5],Ns[1],Ns[2],Ns[3]) end

while ~(is_outer_convergence)

    ρ_in = -ones(Ng) # In-group spectral radius
    if is_CSD
        𝚽E12 = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
    else
        𝚽E12 = Array{Float64}(undef)
    end

    # Loop over energy group
    for ig in range(1,Ng)

        # Calculation of the Legendre components of the source (out-scattering)
        Qlout = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
        if solver_type ∉ [4,5] Qlout = scattering_source(Qlout,𝚽l,Σs[:,:,ig,:],mat,Np,pl,Nm[5],Ns,Ng,ig) end

        # Fixed volumic sources
        Qlout .+= volume_sources[ig,:,:,:,:,:]

        # Calculation of the group flux
        if is_CSD
            if (ig != 1) 𝚽E12 = 𝚽E12 .* ΔE[ig]/ΔE[ig-1] end
            Eg = E[ig]
            ΔEg = ΔE[ig]
            Sg⁻ = S⁻[ig,:]/ΔEg
            Sg⁺ = S⁺[ig,:]/ΔEg
            Sg = S[ig,:,:]/ΔEg
            if solver_type ∈ [2,4]
                Tg = T[ig,:]
            else
                Tg = Vector{Float64}()
                ℳ = Array{Float64}(undef)
            end
        else
            Eg = 0.0
            ΔEg = 0.0
            Sg⁻ = Vector{Float64}()
            Sg⁺ = Vector{Float64}()
            Sg = Vector{Float64}()
            Tg = Vector{Float64}()
            ℳ = Array{Float64}(undef)
        end
        𝚽l[ig,:,:,:,:,:],𝚽E12,ρ_in[ig],Ntot = gn_one_speed(𝚽l[ig,:,:,:,:,:],Qlout,Σtot[ig,:],Σs[:,ig,ig,:],mat,Ndims,ig,Ns,Δs,Np,Nq,pl,pm,Np_surf,𝒪,Nm,isFC,𝒞,ω,I_max,ϵ_max,surface_sources[ig,:,:],is_CSD,solver_type,𝚽E12,Sg⁻,Sg⁺,Sg,Tg,ℳ,𝒜,Ntot,𝒲,Mll,is_SPH,𝒩,boundary_conditions,Np_source,Nv,Mll_surf,Rpq,tiling,gmres_restart,anderson_depth,fold)
    end

    # Verification of convergence in all energy groups
    if is_outer_iteration
        ϵ_out = norm(𝚽l .- 𝚽l⁻) / max(norm(𝚽l), 1e-16)
        𝚽l⁻ = 𝚽l
    end
    if (ϵ_out < ϵ_max || i_out >= I_max) || ~is_outer_iteration
        is_outer_convergence = true
        # Calculate the flux at the cutoff energy
        if is_CSD
            for p in range(1,Np), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[4])
                𝚽cutoff[p,is,ix,iy,iz] = 𝚽E12[p,is,ix,iy,iz]
            end 
        end
    else
        i_out += 1
    end
    
end

# Save flux
flux = Flux_Per_Particle(part)
flux.add_flux(𝚽l)
if is_CSD flux.add_flux_cutoff(𝚽cutoff) end
flux.add_spectral_radius(ρ_in)

return flux

end
