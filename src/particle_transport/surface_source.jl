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
    if get_id(particle) ∉ get_id.(cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Ndims = geometry.get_dimension()
    geometry_type = geometry.get_type()
    norm = 0.0

    surface = source.location

    if solver isa Discrete_Ordinates

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

    elseif solver isa Spherical_Harmonics

        L = min(source.get_legendre_order(),solver.get_legendre_order())
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        Qdims = Ndims
        if polynomial_basis == "legendre"
            Np,_ = half_to_full_range_matrix_legendre(L)
            pl = collect(0:L)
            pm = zeros(Int64,Np)
        elseif polynomial_basis == "spherical-harmonics"
            if Ndims == 1 Qdims = 2 end
            Np,_ = half_to_full_range_matrix_spherical_harmonics(L)
            pl,pm = spherical_harmonics_indices(L)
            error()
        else
            error("Unknown polynomial basis.")
        end
        Lmax = maximum(pl)
    elseif solver isa Galerkin
        L = min(source.get_legendre_order(),solver.get_legendre_order())
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        if polynomial_basis == "legendre"
            error()
        else
            Qdims = 3
            Np,_ = half_to_full_range_matrix_spherical_harmonics(L)
            pl,pm = spherical_harmonics_indices(L)
            Lmax = maximum(pl)
        end
    else
        error("Surface sources are only available with Discrete Ordinates and Spherical Harmonics methods.")
    end

    #----
    # Cartesian 1D geometry
    #----
    if geometry_type == "cartesian" 
        if Ndims == 1

            #----
            # Shifted Legendre polynomials
            #----
            if Qdims == 1

                # Extract source informations
                intensity = source.intensity
                Ωs = source.direction
                μs = Ωs[1]

                # Matrix Initialization
                Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2)
                for ig in range(1,Ng), p in range(1,Np), i in range(1,2)
                    Q[ig,p,i] = 0.0
                end

                # Calculation of the source moments
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

            #----
            # Half-range spherical harmonics
            #----
            else

                # Extract source informations
                intensity = source.intensity
                Ωs = source.direction
                μs = Ωs[1]
                ηs = Ωs[2]
                ξs = Ωs[3]
                ϕs = atan(ξs,ηs)

                # Matrix Initialization
                Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2)
                for ig in range(1,Ng), i in range(1,2), p in range(1,Np)
                    Q[ig,p,i] = 0.0
                end

                # Calculation of the source moments
                norm = intensity
                if solver isa Discrete_Ordinates
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
            Ωs = source.direction
            μs = Ωs[1]
            ηs = Ωs[2]
            ξs = Ωs[3]
            ϕs = atan(ξs,ηs)

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

            # Calculation of the source moments
            norm = 0.0
            if solver isa Discrete_Ordinates
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
                    if surface ∈ ["X-","X+"]
                        for iy in range(1,Ny)
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if surface == "X-" && μs > 0
                                Q[ig,p,1][iy] = intensity * ψlms[l+1][l+m+1] * Δy[iy]
                            elseif surface == "X+" && μs < 0
                                Q[ig,p,2][iy] = intensity * ψlms[l+1][l+m+1] * Δy[iy]
                            end
                            norm += intensity * Δy[iy]
                        end
                    end
                    if surface ∈ ["Y-","Y+"]
                        for ix in range(1,Nx)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if surface == "Y-" && ηs > 0
                                Q[ig,p,3][ix] = intensity * ψlms[l+1][l+m+1] * Δx[ix]
                            elseif surface == "Y+" && ηs < 0
                                Q[ig,p,4][ix] = intensity * ψlms[l+1][l+m+1] * Δx[ix]
                            end
                            norm += intensity * Δx[ix]
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
            Ωs = source.direction
            μs = Ωs[1]
            ηs = Ωs[2]
            ξs = Ωs[3]
            ϕs = atan(ξs,ηs)

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

            # Calculation of the source moments
            norm = 0.0
            if solver isa Discrete_Ordinates
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
                    if surface ∈ ["X-","X+"]
                        for iy in range(1,Ny), iz in range(1,Nz)
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if z[iz] < zmin || z[iz] > zmax continue end
                            if surface == "X-" && μs > 0
                                Q[ig,p,1][iy,iz] = intensity * ψlms[l+1][l+m+1] * Δy[iy] * Δz[iz]
                            elseif surface == "X+" && μs < 0
                                Q[ig,p,2][iy,iz] = intensity * ψlms[l+1][l+m+1] * Δy[iy] * Δz[iz]
                            end
                            norm += intensity * Δy[iy] * Δz[iz]
                        end
                    end
                    if surface ∈ ["Y-","Y+"]
                        for ix in range(1,Nx), iz in range(1,Nz)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if z[iz] < zmin || z[iz] > zmax continue end
                            if surface == "Y-" && ηs > 0
                                Q[ig,p,3][ix,iz] = intensity * ψlms[l+1][l+m+1] * Δx[ix] * Δz[iz]
                            elseif surface == "Y+" && ηs < 0
                                Q[ig,p,4][ix,iz] = intensity * ψlms[l+1][l+m+1] * Δx[ix] * Δz[iz]
                            end
                            norm += intensity * Δx[ix] * Δz[iz]
                        end
                    end
                    if surface ∈ ["Z-","Z+"]
                        for ix in range(1,Nx), iy in range(1,Ny)
                            if x[ix] < xmin || x[ix] > xmax continue end
                            if y[iy] < ymin || y[iy] > ymax continue end
                            if surface == "Z-" && ξs > 0
                                Q[ig,p,5][ix,iy] = intensity * ψlms[l+1][l+m+1] * Δx[ix] * Δy[iy]
                            elseif surface == "Z+" && ξs < 0
                                Q[ig,p,6][ix,iy] = intensity * ψlms[l+1][l+m+1] * Δx[ix] * Δy[iy]
                            end
                            norm += intensity * Δx[ix] * Δy[iy]
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
    error()
    μ⁻,μ⁺,ϕ⁻,ϕ⁺,sb = cartesian_surface_source(b,s)
    #if (b == 3 && (3*π/2 ≤ ϕ ≤ 2*π) && s == -1) || (b == 4 && (3*π/2 ≤ ϕ ≤ 2*π) && s == 1) ϕ -= 2*π end # Adjust for the discontinuity in the azimuthal angle for the Y- boundary
    if μ < μ⁻ || μ > μ⁺ error("Direction cosine is out of bounds for the specified boundary.") end
    if ϕ < ϕ⁻ || ϕ > ϕ⁺ error("Azimuthal angle is out of bounds for the specified boundary.") end
    R_blm = sqrt(2*π/((μ⁺-μ⁻)*(ϕ⁺-ϕ⁻))) * real_half_range_spherical_harmonics_up_to_L(L,-s*sb*((-s*sb-1)/2 + (μ-μ⁻)/(μ⁺-μ⁻)),2*π*(ϕ-ϕ⁻)/(ϕ⁺-ϕ⁻))
    return R_blm
end