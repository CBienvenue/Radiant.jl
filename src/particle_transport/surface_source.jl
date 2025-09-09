"""
    surface_source(particle::Particle,source::Surface_Source,cross_sections::Cross_Sections,
    geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

Prepare the volume source produced by fixed sources for transport calculations.

# Input Argument(s)
- `particle::Particle`: particule type.
- `source::Volume_Source`: volume source information.
- `cross_sections::Cross_Sections`: cross-sections information.
- `geometry::Geometry`: geometry information.
- `discrete_ordinates::Discrete_Ordinates`: discrete_ordinates information.

# Output Argument(s)
- `Q`: boundary sources.
- `norm`: normalization factor.

# Reference(s)
N/A

"""
function surface_source(particle::Particle,source::Surface_Source,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

    #----
    # Initialization
    #----
    if get_id(particle) ∉ get_id.(cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Ndims = geometry.get_dimension()
    geometry_type = geometry.get_type()
    norm = 0.0

    L = source.get_legendre_order()
    quadrature_type = discrete_ordinates.get_quadrature_type()
    N = discrete_ordinates.get_quadrature_order() 
    Qdims = discrete_ordinates.get_quadrature_dimension(Ndims)
    SN_type = discrete_ordinates.get_angular_boltzmann()
    surface = source.location

    # Angular basis
    Ω,w = quadrature(N,quadrature_type,Qdims)
    if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
    Np,Mn,Dn,n⁺_to_n,n_to_n⁺,pl,pm = surface_angular_polynomial_basis(Ω,w,L,SN_type,Qdims,surface)
    Lmax = maximum(pl)

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
                ϕs = atan(Ωs[3],Ωs[2])

                # Matrix Initialization
                Q = Array{Union{Array{Float64},Float64}}(undef,Ng,Np,2)
                for ig in range(1,Ng), i in range(1,2), p in range(1,Np)
                    Q[ig,p,i] = 0.0
                end

                # Calculation of the source moments
                norm = intensity
                ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,μs,ϕs)
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
            if surface == "X-"
                μ⁺ = μs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "X+"
                μ⁺ = -μs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "Y-"
                μ⁺ = ηs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "Y+"
                μ⁺ = -ηs
                ϕ⁺ = atan(ξs,ηs)
            end
            ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,μ⁺,ϕ⁺)
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
            if surface == "X-"
                μ⁺ = μs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "X+"
                μ⁺ = -μs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "Y-"
                μ⁺ = ηs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "Y+"
                μ⁺ = -ηs
                ϕ⁺ = atan(ξs,ηs)
            elseif surface == "Z-"
                μ⁺ = ξs
                ϕ⁺ = atan(ηs,μs)
            elseif surface == "Z+"
                μ⁺ = -ξs
                ϕ⁺ = atan(ηs,μs)
            end
            ψlms = real_half_range_spherical_harmonics_up_to_L(Lmax,μ⁺,ϕ⁺)
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