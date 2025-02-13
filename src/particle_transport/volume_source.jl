"""
    volume_source(particle::Particle,source::Volume_Source,cross_sections::Cross_Sections,
    geometry::Geometry)

Prepare the volume source produced by fixed sources for transport calculations.

# Input Argument(s)
- `particle::Particle`: particule type.
- `source::Volume_Source`: volume source information.
- `cross_sections::Cross_Sections`: cross-sections information.
- `geometry::Geometry`: geometry information.

# Output Argument(s)
- `Q`: source density.
- `norm`: normalization factor.

# Reference(s)
N/A

"""
function volume_source(particle::Particle,source::Volume_Source,cross_sections::Cross_Sections,geometry::Geometry)

    if get_id(particle) ∉ get_id.(cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity

    if geometry.type == "cartesian"

    #----
    # Cartesian 1D geometry
    #----
    if geometry.dimension == 1
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        Q = zeros(Ng,1,1,Nx,1,1)
        norm = 0
        for ig in range(1,Ng)
        if energy_groups[ig] == 1 
            for ix in range(1,Nx)
                if xmin <= x[ix] && xmax >= x[ix]
                    Q[ig,1,1,ix,1,1] += intensity * Δx[ix]
                    norm = norm + intensity * Δx[ix]
                end
            end
        end
        end
    #----
    # Cartesian 2D geometry
    #----
    elseif geometry.dimension == 2
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        y = geometry.voxels_position["y"]
        Δy = geometry.voxels_width["y"]
        Ny = geometry.number_of_voxels["y"]
        ymin = source.boundaries["y"][1]
        ymax = source.boundaries["y"][2]
        Q = zeros(Ng,1,1,Nx,Ny,1)
        norm = 0
        for ig in range(1,Ng)
        if energy_groups[ig] == 1 
            for ix in range(1,Nx), iy in range(1,Ny)
                if xmin <= x[ix] && xmax >= x[ix] && ymin <= y[iy] && ymax >= y[iy] 
                    Q[ig,1,1,ix,iy,1] += intensity * Δx[ix] * Δy[iy]
                    norm = norm + intensity * Δx[ix] * Δy[iy]
                end
            end
        end
        end
    #----
    # Cartesian 3D geometry
    #----
    elseif geometry.dimension == 3
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        y = geometry.voxels_position["y"]
        Δy = geometry.voxels_width["y"]
        Ny = geometry.number_of_voxels["y"]
        ymin = source.boundaries["y"][1]
        ymax = source.boundaries["y"][2]
        z = geometry.voxels_position["z"]
        Δz = geometry.voxels_width["z"]
        Nz = geometry.number_of_voxels["z"]
        zmin = source.boundaries["z"][1]
        zmax = source.boundaries["z"][2]
        Q = zeros(Ng,1,1,Nx,Ny,Nz)
        norm = 0
        for ig in range(1,Ng)
        if energy_groups[ig] == 1 
            for ix in range(1,Nx), iy in range(1,Ny), iz in range(1,Nz)
                if xmin <= x[ix] && xmax >= x[ix] && ymin <= y[iy] && ymax >= y[iy] && zmin <= z[iz] && zmax >= z[iz] 
                    Q[ig,1,1,ix,iy,iz] += intensity * Δx[ix] * Δy[iy] * Δz[iz]
                    norm = norm + intensity * Δx[ix] * Δy[iy] * Δz[iz]
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