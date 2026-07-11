"""
    surface_source(particle::Particle,source::Surface_Source,cross_sections::Cross_Sections,
    geometry::Geometry,solver::Solver)

Prepare the volume source produced by fixed sources for transport calculations.

# Input Argument(s)
- `particle::Particle`: particule type.
- `source::Volume_Source`: volume source information.
- `cross_sections::Cross_Sections`: cross-sections information.
- `geometry::Geometry`: geometry information.
- `solver::Solver`: solver information.

# Output Argument(s)
- `Q`: boundary sources.
- `norm`: normalization factor.

# Reference(s)
N/A

"""
function surface_source(particle::Particle,source::Surface_Source,cross_sections::Cross_Sections,geometry::Geometry,solver::Solver)

    #----
    # Initialization
    #----
    if get_tag(particle) ∉ get_tag.(cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    index = findfirst(x -> get_tag(x) == get_tag(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Ndims = geometry.get_dimension()
    geometry_type = geometry.get_type()
    norm = 0.0
    surface = uppercase(source.location)
    if ~ismissing(source.angular_moments) && Ndims != 1 && surface ∉ ["X-","X+"]
        error("Distributed surface sources (set_angular_moments) are only available on the X- and X+ faces in 2D/3D.")
    end

    if solver isa SN

        L = min(source.get_legendre_order(),solver.get_legendre_order())
        quadrature_type = solver.get_quadrature_type()
        N = solver.get_quadrature_order() 
        Qdims = solver.get_quadrature_dimension(Ndims)
        SN_type = solver.get_angular_boltzmann()

        # Angular basis
        Ω,w = quadrature(N,quadrature_type,Qdims)
        if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
        Np,_,_,_,_,pl,pm = surface_angular_polynomial_basis(Ω,w,L,SN_type,Qdims,surface)
        Lmax = maximum(pl)

    elseif solver isa GN
        L = min(source.get_legendre_order(),solver.get_legendre_order())
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        if polynomial_basis == "legendre"
            Qdims = 1
            Np,_ = half_to_full_range_matrix_legendre(L)
            pl = collect(0:L)
            pm = zeros(Int64,Np)
            Lmax = maximum(pl)
        else
            Qdims = 3
            Np,_ = half_to_full_range_matrix_spherical_harmonics(L)
            pl,pm = spherical_harmonics_indices(L)
            Lmax = maximum(pl)
        end
    elseif solver isa CP
        # Azimuthally-symmetric (m = 0) half-range Legendre surface expansion, truncated at the
        # solver's surface order Nν. The incoming boundary flux is represented by its moments in
        # the same basis R̄_ℓ(μ̂) = √(2ℓ+1) P_ℓ(2μ̂-1) used by the CP boundary-coupling matrices.
        if Ndims != 1 error("The CP solver only supports surface sources in 1D Cartesian geometry.") end
        L = min(source.get_legendre_order(),solver.get_surface_order())
        Qdims = 1
        Np = L+1
        pl = collect(0:L)
        pm = zeros(Int64,Np)
        Lmax = L
    else
        error("Surface sources are only available with Discrete Ordinates, Spherical Harmonics and Collision Probability methods.")
    end

    #----
    # Cartesian 1D geometry
    #----
    if geometry_type == "cartesian"
        if Ndims == 1

            is_distributed = ~ismissing(source.angular_moments)
            if is_distributed && surface ∉ ["X-","X+"] error("Distributed surface sources (set_angular_moments) are only available on the X- and X+ faces.") end

            #----
            # Shifted Legendre polynomials
            #----
            if Qdims == 1

                # Extract source informations
                intensity = source.intensity

                # Matrix Initialization
                Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2)
                for ig in range(1,Ng), p in range(1,Np), i in range(1,2)
                    Q[ig,p,i] = 0.0
                end

                # Calculation of the source moments
                if is_distributed
                    # Distributed source: the user-provided half-range moments of the
                    # incident angular flux are already in the surface basis R̄ₚ(μ̂).
                    g_mom = source.angular_moments
                    norm = intensity * g_mom[1]
                    ib = surface == "X-" ? 1 : 2
                    for ig in range(1,Ng), p in range(1,Np)
                        if ig == source.energy_group && pl[p] ≤ length(g_mom)-1
                            Q[ig,p,ib] = intensity * g_mom[pl[p]+1]
                        end
                    end
                else
                    μs = source.direction[1]
                    norm = intensity
                    Pls = half_range_legendre_polynomials_up_to_L(Lmax,abs(μs))
                    for ig in range(1,Ng), p in range(1,Np)
                        if ig == source.energy_group
                            l = pl[p]
                            if surface == "X-"
                                if (μs > 0) Q[ig,p,1] = intensity * Pls[l+1] end
                            elseif surface == "X+"
                                if (μs < 0) Q[ig,p,2] = intensity * Pls[l+1] end
                            end
                        end
                    end
                end

            #----
            # Half-range spherical harmonics
            #----
            else

                # Extract source informations
                intensity = source.intensity

                # Matrix Initialization
                Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2)
                for ig in range(1,Ng), i in range(1,2), p in range(1,Np)
                    Q[ig,p,i] = 0.0
                end

                # Calculation of the source moments
                if is_distributed
                    # Distributed (azimuthally-symmetric) source: the m = 0 half-range
                    # spherical harmonics are R̄ₗ(μ̂)/√(2π) on the x-faces, so the surface
                    # moments are the user-provided half-range Legendre moments ÷ √(2π).
                    # The GN boundary basis (boundary_real_half_range_...) maps the face
                    # box with the REVERSED argument μ̂ → 1-μ̂, so its odd-ℓ members carry
                    # an extra (-1)^ℓ.
                    g_mom = source.angular_moments
                    norm = intensity * g_mom[1]
                    ib = surface == "X-" ? 1 : 2
                    for ig in range(1,Ng), p in range(1,Np)
                        if ig == source.energy_group && pm[p] == 0 && pl[p] ≤ length(g_mom)-1
                            sgn = (solver isa GN && isodd(pl[p])) ? -1.0 : 1.0
                            Q[ig,p,ib] = sgn * intensity * g_mom[pl[p]+1] / sqrt(2*π)
                        end
                    end
                else
                    Ωs = source.direction
                    μs = Ωs[1]
                    ηs = Ωs[2]
                    ξs = Ωs[3]
                    ϕs = atan(ξs,ηs)
                    norm = intensity
                    if solver isa SN
                        ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,abs(μs),ϕs)
                    else
                        b = cartesian_boundary_index(surface)
                        ψlms = boundary_real_half_range_spherical_harmonics_up_to_L(L,b,-1,μs,ϕs)
                    end
                    for ig in range(1,Ng)
                        if ~(ig == source.energy_group) continue end
                        for p in range(1,Np)
                            l = pl[p]
                            m = pm[p]
                            if surface == "X-"
                                if (μs > 0) Q[ig,p,1] = intensity * ψlms[l+1][l+m+1] end
                            elseif surface == "X+"
                                if (μs < 0) Q[ig,p,2] = intensity * ψlms[l+1][l+m+1] end
                            end
                        end
                    end
                end
            end

        elseif Ndims == 2

            # Geometry information
            x = geometry.voxels_position["x"]
            Δx = geometry.voxels_width["x"]
            Nx = geometry.number_of_voxels["x"]
            y = geometry.voxels_position["y"]
            Δy = geometry.voxels_width["y"]
            Ny = geometry.number_of_voxels["y"]
            if surface ∈ ["X-","X+"]
                ymin = source.boundaries["y"][1]
                ymax = source.boundaries["y"][2]
            elseif surface ∈ ["Y-","Y+"]
                xmin = source.boundaries["x"][1]
                xmax = source.boundaries["x"][2]
            end

            # Extract source informations
            intensity = source.intensity
            is_distributed = ~ismissing(source.angular_moments)
            if ~is_distributed
                Ωs = source.direction
                μs = Ωs[1]
                ηs = Ωs[2]
                ξs = Ωs[3]
                ϕs = atan(ξs,ηs)
            end

            # Matrix Initialization
            Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,4)
            for ig in range(1,Ng), p in range(1,Np)
                for i in range(1,2)
                    Q[ig,p,i] = zeros(Ny)
                end
                for i in range(3,4)
                    Q[ig,p,i] = zeros(Nx)
                end
            end

            # Angular value per moment slot: half-range harmonics at the δ direction,
            # or the m = 0 collapse of the user moments (÷√(2π), x-faces span the full
            # azimuth) for a distributed source.
            aval = zeros(Np)
            if is_distributed
                # The GN boundary basis (boundary_real_half_range_...) maps the face
                # box with the REVERSED argument μ̂ → 1-μ̂ on the x-faces, so its odd-ℓ
                # members carry an extra (-1)^ℓ; the SN basis is unreversed.
                g_mom = source.angular_moments
                for p in range(1,Np)
                    if pm[p] == 0 && pl[p] ≤ length(g_mom)-1
                        sgn = (solver isa GN && isodd(pl[p])) ? -1.0 : 1.0
                        aval[p] = sgn*g_mom[pl[p]+1]/sqrt(2*π)
                    end
                end
            else
                if solver isa SN
                    ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,abs(μs),ϕs)
                else
                    b = cartesian_boundary_index(surface)
                    ψlms = boundary_real_half_range_spherical_harmonics_up_to_L(L,b,-1,μs,ϕs)
                end
                for p in range(1,Np)
                    aval[p] = ψlms[pl[p]+1][pl[p]+pm[p]+1]
                end
            end
            g0fac = is_distributed ? source.angular_moments[1] : 1.0
            in_xm = is_distributed || μs > 0
            in_xp = is_distributed || μs < 0

            # Calculation of the source moments
            norm = 0.0
            for ig in range(1,Ng)
                if ~(ig == source.energy_group) continue end
                for p in range(1,Np)
                    # The boundary moments are angular-flux VALUES on each face cell
                    # (the sweeps consume them as incoming cell-edge fluxes); the
                    # transverse widths enter only the source normalization.
                    if surface ∈ ["X-","X+"]
                        for iy in range(1,Ny)
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if surface == "X-" && in_xm
                                Q[ig,p,1][iy] = intensity * aval[p]
                            elseif surface == "X+" && in_xp
                                Q[ig,p,2][iy] = intensity * aval[p]
                            end
                            if p == 1 norm += intensity * g0fac * Δy[iy] end
                        end
                    end
                    if surface ∈ ["Y-","Y+"]
                        for ix in range(1,Nx)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if surface == "Y-" && ηs > 0
                                Q[ig,p,3][ix] = intensity * aval[p]
                            elseif surface == "Y+" && ηs < 0
                                Q[ig,p,4][ix] = intensity * aval[p]
                            end
                            if p == 1 norm += intensity * Δx[ix] end
                        end
                    end
                end
            end

        elseif Ndims == 3
            
            # Geometry information
            x = geometry.voxels_position["x"]
            Δx = geometry.voxels_width["x"]
            Nx = geometry.number_of_voxels["x"]
            y = geometry.voxels_position["y"]
            Δy = geometry.voxels_width["y"]
            Ny = geometry.number_of_voxels["y"]
            z = geometry.voxels_position["z"]
            Δz = geometry.voxels_width["z"]
            Nz = geometry.number_of_voxels["z"]
            if surface ∈ ["X-","X+"]
                ymin = source.boundaries["y"][1]
                ymax = source.boundaries["y"][2]
                zmin = source.boundaries["z"][1]
                zmax = source.boundaries["z"][2]
            elseif surface ∈ ["Y-","Y+"]
                xmin = source.boundaries["x"][1]
                xmax = source.boundaries["x"][2]
                zmin = source.boundaries["z"][1]
                zmax = source.boundaries["z"][2]
            elseif surface ∈ ["Z-","Z+"]
                xmin = source.boundaries["x"][1]
                xmax = source.boundaries["x"][2]
                ymin = source.boundaries["y"][1]
                ymax = source.boundaries["y"][2]
            end

            # Extract source informations
            intensity = source.intensity
            is_distributed = ~ismissing(source.angular_moments)
            if ~is_distributed
                Ωs = source.direction
                μs = Ωs[1]
                ηs = Ωs[2]
                ξs = Ωs[3]
                ϕs = atan(ξs,ηs)
            end

            # Matrix Initialization
            Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,6)
            for ig in range(1,Ng), p in range(1,Np)
                for i in range(1,2)
                    Q[ig,p,i] = zeros(Ny,Nz)
                end
                for i in range(3,4)
                    Q[ig,p,i] = zeros(Nx,Nz)
                end
                for i in range(5,6)
                    Q[ig,p,i] = zeros(Nx,Ny)
                end
            end

            # Angular value per moment slot: half-range harmonics at the δ direction,
            # or the m = 0 collapse of the user moments (÷√(2π), x-faces span the full
            # azimuth) for a distributed source.
            aval = zeros(Np)
            if is_distributed
                # The GN boundary basis (boundary_real_half_range_...) maps the face
                # box with the REVERSED argument μ̂ → 1-μ̂ on the x-faces, so its odd-ℓ
                # members carry an extra (-1)^ℓ; the SN basis is unreversed.
                g_mom = source.angular_moments
                for p in range(1,Np)
                    if pm[p] == 0 && pl[p] ≤ length(g_mom)-1
                        sgn = (solver isa GN && isodd(pl[p])) ? -1.0 : 1.0
                        aval[p] = sgn*g_mom[pl[p]+1]/sqrt(2*π)
                    end
                end
            else
                if solver isa SN
                    ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,abs(μs),ϕs)
                else
                    b = cartesian_boundary_index(surface)
                    ψlms = boundary_real_half_range_spherical_harmonics_up_to_L(L,b,-1,μs,ϕs)
                end
                for p in range(1,Np)
                    aval[p] = ψlms[pl[p]+1][pl[p]+pm[p]+1]
                end
            end
            g0fac = is_distributed ? source.angular_moments[1] : 1.0
            in_xm = is_distributed || μs > 0
            in_xp = is_distributed || μs < 0

            # Calculation of the source moments
            norm = 0.0
            for ig in range(1,Ng)
                if ~(ig == source.energy_group) continue end
                for p in range(1,Np)
                    # The boundary moments are angular-flux VALUES on each face cell
                    # (the sweeps consume them as incoming cell-edge fluxes); the
                    # transverse widths enter only the source normalization.
                    if surface ∈ ["X-","X+"]
                        for iy in range(1,Ny), iz in range(1,Nz)
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if z[iz] < zmin || z[iz] > zmax continue end
                            if surface == "X-" && in_xm
                                Q[ig,p,1][iy,iz] = intensity * aval[p]
                            elseif surface == "X+" && in_xp
                                Q[ig,p,2][iy,iz] = intensity * aval[p]
                            end
                            if p == 1 norm += intensity * g0fac * Δy[iy] * Δz[iz] end
                        end
                    end
                    if surface ∈ ["Y-","Y+"]
                        for ix in range(1,Nx), iz in range(1,Nz)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if z[iz] < zmin || z[iz] > zmax continue end
                            if surface == "Y-" && ηs > 0
                                Q[ig,p,3][ix,iz] = intensity * aval[p]
                            elseif surface == "Y+" && ηs < 0
                                Q[ig,p,4][ix,iz] = intensity * aval[p]
                            end
                            if p == 1 norm += intensity * Δx[ix] * Δz[iz] end
                        end
                    end
                    if surface ∈ ["Z-","Z+"]
                        for ix in range(1,Nx), iy in range(1,Ny)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if surface == "Z-" && ξs > 0
                                Q[ig,p,5][ix,iy] = intensity * aval[p]
                            elseif surface == "Z+" && ξs < 0
                                Q[ig,p,6][ix,iy] = intensity * aval[p]
                            end
                            if p == 1 norm += intensity * Δx[ix] * Δy[iy] end
                        end
                    end
                end
            end

        else
            error("Cartesian geometries are available only in 1D, 2D and 3D geometries.")
        end
    else
        error("Surface sources are only available in Cartesian geometries.")
    end
    return Q, norm
end

function cartesian_boundary_index(location::String)
    if location == "X-"
        return 1
    elseif location == "X+"
        return 2
    elseif location == "Y-"
        return 3
    elseif location == "Y+"
        return 4
    elseif location == "Z-"
        return 5
    elseif location == "Z+"
        return 6
    else
        error("Unknown surface source location.")
    end
end

function cartesian_surface_source(b::Int64,s::Int64)
    if b < 1 || b > 6 error("Boundary index must be between 1 and 6.") end
    if b == 1
        μ⁻ = -(s+1)/2
        μ⁺ = (-s+1)/2
        ϕ⁻ = 0
        ϕ⁺ = 2π
        sb = -1
    elseif b == 2
        μ⁻ = (s-1)/2
        μ⁺ = (s+1)/2
        ϕ⁻ = 0
        ϕ⁺ = 2π
        sb = 1
    elseif b == 3
        μ⁻ = -1
        μ⁺ = 1
        ϕ⁻ = π/2 + (s-1)/2 * π
        ϕ⁺ = 3*π/2 + (s-1)/2 * π
        sb = -1
    elseif b == 4
        μ⁻ = -1
        μ⁺ = 1
        ϕ⁻ = -π/2 - (s-1)/2 * π
        ϕ⁺ = π/2 - (s-1)/2 * π
        sb = 1
    elseif b == 5
        μ⁻ = -1
        μ⁺ = 1
        ϕ⁻ = π + (s-1)/2 * π
        ϕ⁺ = 2*π + (s-1)/2 * π
        sb = -1
    elseif b == 6
        μ⁻ = -1
        μ⁺ = 1
        ϕ⁻ = 0 - (s-1)/2 * π
        ϕ⁺ = π - (s-1)/2 * π
        sb = 1
    end
    return μ⁻, μ⁺, ϕ⁻, ϕ⁺, sb
end

function boundary_real_half_range_spherical_harmonics_up_to_L(L::Int64,b::Int64,s::Int64,μ::Float64,ϕ::Float64)
    μ⁻,μ⁺,ϕ⁻,ϕ⁺,sb = cartesian_surface_source(b,s)
    #if (b == 3 && (3*π/2 ≤ ϕ ≤ 2*π) && s == -1) || (b == 4 && (3*π/2 ≤ ϕ ≤ 2*π) && s == 1) ϕ -= 2*π end # Adjust for the discontinuity in the azimuthal angle for the Y- boundary
    if μ < μ⁻ || μ > μ⁺ error("Direction cosine is out of bounds for the specified boundary.") end
    if ϕ < ϕ⁻ || ϕ > ϕ⁺ error("Azimuthal angle is out of bounds for the specified boundary.") end
    R_blm = sqrt(2*π/((μ⁺-μ⁻)*(ϕ⁺-ϕ⁻))) * real_half_range_spherical_harmonics_up_to_L(L,-s*sb*((-s*sb-1)/2 + (μ-μ⁻)/(μ⁺-μ⁻)),2*π*(ϕ-ϕ⁻)/(ϕ⁺-ϕ⁻))
    return R_blm
end