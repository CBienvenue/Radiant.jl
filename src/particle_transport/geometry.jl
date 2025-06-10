"""
    geometry(geo::Geometry,cs::Cross_Sections)

Compute the geometry properties for transport calculations. 

# Input Argument(s)
- `geo::Geometry`: geometry informations.
- `cs::Cross_Sections`: cross section informations.

# Output Argument(s)
- `geo::Geometry`: geometry informations.

# Reference(s)
N/A

"""
function geometry(geo::Geometry,cs::Cross_Sections)

    #----
    # 1D Cartesian geometry
    #----

    if geo.type == "cartesian" && geo.dimension == 1

    # Extracting data
    Nrx = geo.number_of_regions["x"]
    xrb = geo.region_boundaries["x"]
    Nrx_voxels = geo.voxels_per_region["x"]
    idr = geo.material_per_region
    Nmat = length(cs.materials)
    id_mat = [mat.get_id() for mat in cs.materials]

    # Compute positions and widths along x-axis
    Nx = sum(Nrx_voxels)
    xb = zeros(Nx+1)
    Δx = zeros(Nx)
    x = zeros(Nx)
    ix = 0
    xp = xrb[1]
    for ir in range(1,Nrx)
        Δxi = (xrb[ir+1]-xrb[ir])/Nrx_voxels[ir]
        for iv in range(1,Nrx_voxels[ir])
            xp += Δxi
            ix += 1
            xb[ix+1] = xp
            x[ix] = xp-Δxi/2
            Δx[ix] = Δxi
        end
    end

    # Compute volumes
    vol = Δx

    # Set material matrix
    mat=zeros(Int64,Nx,1,1)
    ix = 0
    for ir in range(1,Nrx)
        for iv in range(1,Nrx_voxels[ir])
            ix += 1
            index = findfirst(x->x==idr[ir].get_id(),id_mat)
            if isnothing(index) error("No cross section information for material $(idr[ir].get_id()).") end
            mat[ix] = index
        end
    end

    # Save data
    geo.voxels_width["x"] = Δx
    geo.voxels_position["x"] = x
    geo.voxels_boundaries["x"] = xb
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    #----
    # 2D Cartesian geometry
    #----

    elseif geo.type == "cartesian" && geo.dimension == 2

    # Extracting data
    Nrx = geo.number_of_regions["x"]
    Nry = geo.number_of_regions["y"]
    xrb = geo.region_boundaries["x"]
    yrb = geo.region_boundaries["y"]
    Nrx_voxels = geo.voxels_per_region["x"]
    Nry_voxels = geo.voxels_per_region["y"]
    idr = geo.material_per_region
    Nmat = length(cs.materials)
    id_mat = [mat.get_id() for mat in cs.materials]

    # Compute positions and widths along x-axis
    Nx = sum(Nrx_voxels)
    xb = zeros(Nx+1)
    Δx = zeros(Nx)
    x = zeros(Nx)
    irx = zeros(Int64,Nx)
    ix = 0
    xp = xrb[1]
    for ir in range(1,Nrx)
        Δxi = (xrb[ir+1]-xrb[ir])/Nrx_voxels[ir]
        for iv in range(1,Nrx_voxels[ir])
            xp += Δxi
            ix += 1
            xb[ix+1] = xp
            x[ix] = xp-Δxi/2
            Δx[ix] = Δxi
            irx[ix] = ir
        end
    end

    # Compute positions and widths along y-axis
    Ny = sum(Nry_voxels)
    yb = zeros(Ny+1)
    Δy = zeros(Ny)
    y = zeros(Ny)
    iry = zeros(Int64,Ny)
    iy = 0
    yp = yrb[1]
    for ir in range(1,Nry)
        Δyi = (yrb[ir+1]-yrb[ir])/Nry_voxels[ir]
        for iv in range(1,Nry_voxels[ir])
            yp += Δyi
            iy += 1
            yb[iy+1] = yp
            y[iy] = yp-Δyi/2
            Δy[iy] = Δyi
            iry[iy] = ir
        end
    end

    # Compute volumes
    vol = zeros(Nx,Ny)
    for ix in range(1,Nx), iy in range(1,Ny)
        vol[ix,iy] = Δx[ix] * Δy[iy]
    end

    # Set material matrix
    mat=zeros(Int64,Nx,Ny,1)
    for ix in range(1,Nx), iy in range(1,Ny)
        index = findfirst(x->x==idr[irx[ix],iry[iy]].get_id(),id_mat)
        if isnothing(index) error("No cross section information for material $(idr[irx[ix],iry[iy]].get_id()).") end
        mat[ix,iy] = index
    end

    # Save data
    geo.voxels_width["x"] = Δx
    geo.voxels_width["y"] = Δy
    geo.voxels_position["x"] = x
    geo.voxels_position["y"] = y
    geo.voxels_boundaries["x"] = xb
    geo.voxels_boundaries["y"] = yb
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    #----
    # 3D Cartesian geometry
    #----

    elseif geo.type == "cartesian" && geo.dimension == 3

    # Extracting data
    Nrx = geo.number_of_regions["x"]
    Nry = geo.number_of_regions["y"]
    Nrz = geo.number_of_regions["z"]
    xrb = geo.region_boundaries["x"]
    yrb = geo.region_boundaries["y"]
    zrb = geo.region_boundaries["z"]
    Nrx_voxels = geo.voxels_per_region["x"]
    Nry_voxels = geo.voxels_per_region["y"]
    Nrz_voxels = geo.voxels_per_region["z"]
    idr = geo.material_per_region
    Nmat = length(cs.materials)
    id_mat = [mat.get_id() for mat in cs.materials]

    # Compute positions and widths along x-axis
    Nx = sum(Nrx_voxels)
    xb = zeros(Nx+1)
    Δx = zeros(Nx)
    x = zeros(Nx)
    irx = zeros(Int64,Nx)
    ix = 0
    xp = xrb[1]
    for ir in range(1,Nrx)
        Δxi = (xrb[ir+1]-xrb[ir])/Nrx_voxels[ir]
        for iv in range(1,Nrx_voxels[ir])
            xp += Δxi
            ix += 1
            xb[ix+1] = xp
            x[ix] = xp-Δxi/2
            Δx[ix] = Δxi
            irx[ix] = ir
        end
    end

    # Compute positions and widths along y-axis
    Ny = sum(Nry_voxels)
    yb = zeros(Ny+1)
    Δy = zeros(Ny)
    y = zeros(Ny)
    iry = zeros(Int64,Ny)
    iy = 0
    yp = yrb[1]
    for ir in range(1,Nry)
        Δyi = (yrb[ir+1]-yrb[ir])/Nry_voxels[ir]
        for iv in range(1,Nry_voxels[ir])
            yp += Δyi
            iy += 1
            yb[iy+1] = yp
            y[iy] = yp-Δyi/2
            Δy[iy] = Δyi
            iry[iy] = ir
        end
    end

    # Compute positions and widths along z-axis
    Nz = sum(Nrz_voxels)
    zb = zeros(Nz+1)
    Δz = zeros(Nz)
    z = zeros(Nz)
    irz = zeros(Int64,Nz)
    iz = 0
    zp = zrb[1]
    for ir in range(1,Nrz)
        Δzi = (zrb[ir+1]-zrb[ir])/Nrz_voxels[ir]
        for iv in range(1,Nrz_voxels[ir])
            zp += Δzi
            iz += 1
            zb[iz+1] = zp
            z[iz] = zp-Δzi/2
            Δz[iz] = Δzi
            irz[iz] = ir
        end
    end

    # Compute volumes
    vol = zeros(Nx,Ny,Nz)
    for ix in range(1,Nx), iy in range(1,Ny), iz in range(1,Nz)
        vol[ix,iy,iz] = Δx[ix] * Δy[iy] * Δz[iz]
    end

    # Set material matrix
    mat=zeros(Int64,Nx,Ny,Nz)
    for ix in range(1,Nx), iy in range(1,Ny), iz in range(1,Nz)
        index = findfirst(x->x==idr[irx[ix],iry[iy],irz[iz]].get_id(),id_mat)
        if isnothing(index) error("No cross section information for material $(idr[irx[ix],iry[iy],irz[iz]].get_id()).") end
        mat[ix,iy,iz] = index
    end

    # Save data
    geo.voxels_width["x"] = Δx
    geo.voxels_width["y"] = Δy
    geo.voxels_width["z"] = Δz
    geo.voxels_position["x"] = x
    geo.voxels_position["y"] = y
    geo.voxels_position["z"] = z
    geo.voxels_boundaries["x"] = xb
    geo.voxels_boundaries["y"] = yb
    geo.voxels_boundaries["z"] = zb
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    else
        error("Undefined type of geometry.")
    end
    return geo
end