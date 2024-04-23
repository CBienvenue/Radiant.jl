"""
    geometry(geo::Geometry,cs::Cross_Sections)

Compute the geometry properties for transport calculations. 

See also [`transport`](@ref).

# Input Argument(s)
- 'geo::Geometry': geometry informations.
- 'cs::Cross_Sections': cross section informations.

# Output Argument(s)
- 'geo::Geometry': geometry informations.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function geometry(geo::Geometry,cs::Cross_Sections)

    #----
    # 1D Cartesian geometry
    #----

    if geo.type == "cartesian" && geo.dimension == 1

    # Extracting data

    Nrx           = geo.number_of_regions["x"]
    xB            = geo.region_boundaries["x"]
    xVoxels       = geo.voxels_per_region["x"]
    idMat         = geo.material_per_region

    Nmat = length(cs.materials)
    materialNames = Vector{String}()
    for n in range(1,Nmat)
        push!(materialNames,cs.materials[n].name)
    end

    # Calculations of the volume and material array

    nVoxels = sum(xVoxels)
    vol=zeros(nVoxels)
    mat=zeros(Int64,nVoxels,1,1)

    i=1
    while i <= Nrx
        if i==1
            rx = range(1,xVoxels[i])
        else
            rx = range(sum(xVoxels[1:i-1])+1,sum(xVoxels[1:i-1])+xVoxels[i])
        end
        for ix in rx
            vol[ix] = (xB[i+1]-xB[i])/xVoxels[i]
            index = findfirst(x->x==idMat[i],materialNames)
            if isnothing(index) error("No cross section information for material ",idMat[i],".") end
            mat[ix] = index
        end
        i += 1
    end

    # Calculations of the x-coordinates of each voxels

    k=1
    xp = zeros(Int64,nVoxels)
    x = zeros(nVoxels)

    for i in range(1,Nrx)
        for j in range(1,xVoxels[i])
            xp[k]=i
            k=k+1
        end
    end

    for ix in range(1,nVoxels)
        stepx = (xB[xp[ix]+1]-xB[xp[ix]])/xVoxels[xp[ix]]
        if xp[ix] == 1
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-1) * stepx
        else
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-sum(xVoxels[1:xp[ix]-1])-1) * stepx
        end
    end

    # Save data

    geo.voxels_width["x"] = vol
    geo.voxels_position["x"] = x
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    #----
    # 2D Cartesian geometry
    #----

    elseif geo.type == "cartesian" && geo.dimension == 2

    # Extracting data

    Nrx           = geo.number_of_regions["x"]
    xB            = geo.region_boundaries["x"]
    xVoxels       = geo.voxels_per_region["x"]
    Nry           = geo.number_of_regions["y"]
    yB            = geo.region_boundaries["y"]
    yVoxels       = geo.voxels_per_region["y"]
    idMat         = geo.material_per_region

    Nmat = length(cs.materials)
    materialNames = Vector{String}()
    for n in range(1,Nmat)
        push!(materialNames,cs.materials[n].name)
    end

    # Calculations of the volume and material array

    nVoxels = [sum(xVoxels),sum(yVoxels)]
    vol = zeros(nVoxels[1],nVoxels[2])
    mat = zeros(Int64,nVoxels[1],nVoxels[2],1)
    Δx = zeros(nVoxels[1])
    Δy = zeros(nVoxels[2])

    # Voxels x-width
    i=1
    while i <= Nrx
        if i==1
            rx = range(1,xVoxels[i])
        else
            rx = range(sum(xVoxels[1:i-1])+1,sum(xVoxels[1:i-1])+xVoxels[i])
        end
        for ix in rx
            Δx[ix] = (xB[i+1]-xB[i])/xVoxels[i]
        end
        i += 1
    end

    # Voxels y-width
    j=1
    while j <= Nry
        if j==1
            ry = range(1,yVoxels[j])
        else
            ry = range(sum(yVoxels[1:j-1])+1,sum(yVoxels[1:j-1])+yVoxels[j])
        end
        for iy in ry
            Δy[iy] = (yB[j+1]-yB[j])/yVoxels[j]
        end
        j += 1
    end

    # Volumes
    for ix in range(1,nVoxels[1])
        for iy in range(1,nVoxels[2])
            vol[ix,iy] = Δx[ix] * Δy[iy]
        end
    end

    # Material matrix
    i=1
    j=1
    while i <= Nrx
        if i==1
            rx = range(1,xVoxels[i])
        else
            rx = range(sum(xVoxels[1:i-1])+1,sum(xVoxels[1:i-1])+xVoxels[i])
        end
        while j <= Nry
            if j==1
                ry = range(1,yVoxels[j])
            else
                ry = range(sum(yVoxels[1:j-1])+1,sum(yVoxels[1:j-1])+yVoxels[j])
            end
            for ix in rx
                for iy in ry
                    index = findfirst(x->x==idMat[i,j],materialNames)
                    if isnothing(index) error("No cross section information for material ",idMat[i,j],".") end
                    mat[ix,iy] = index
                end
            end
            j += 1
        end
        j = 1
        i += 1
    end

    # Calculations of the x-coordinates of each voxels

    k=1
    xp = zeros(Int64,nVoxels[1])
    x = zeros(nVoxels[1])

    for i in range(1,Nrx)
        for j in range(1,xVoxels[i])
            xp[k]=i
            k=k+1
        end
    end

    for ix in range(1,nVoxels[1])
        stepx = (xB[xp[ix]+1]-xB[xp[ix]])/xVoxels[xp[ix]]
        if xp[ix] == 1
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-1) * stepx
        else
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-sum(xVoxels[1:xp[ix]-1])-1) * stepx
        end
    end

    # Calculations of the y-coordinates of each voxels

    k=1
    yp = zeros(Int64,nVoxels[2])
    y = zeros(nVoxels[2])

    for i in range(1,Nry)
        for j in range(1,yVoxels[i])
            yp[k]=i
            k=k+1
        end
    end

    for iy in range(1,nVoxels[2])
        stepy = (yB[yp[iy]+1]-yB[yp[iy]])/yVoxels[yp[iy]]
        if yp[iy] == 1
            y[iy] = yB[yp[iy]] + 0.5 * stepy + (iy-1) * stepy
        else
            y[iy] = yB[yp[iy]] + 0.5 * stepy + (iy-sum(yVoxels[1:yp[iy]-1])-1) * stepy
        end
    end

    # Save data

    geo.voxels_width["x"] = Δx
    geo.voxels_width["y"] = Δy
    geo.voxels_position["x"] = x
    geo.voxels_position["y"] = y
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    #----
    # 3D Cartesian geometry
    #----

    elseif geo.type == "cartesian" && geo.dimension == 3

    # Extracting data

    Nrx           = geo.number_of_regions["x"]
    xB            = geo.region_boundaries["x"]
    xVoxels       = geo.voxels_per_region["x"]
    Nry           = geo.number_of_regions["y"]
    yB            = geo.region_boundaries["y"]
    yVoxels       = geo.voxels_per_region["y"]
    Nrz           = geo.number_of_regions["z"]
    zB            = geo.region_boundaries["z"]
    zVoxels       = geo.voxels_per_region["z"]
    idMat         = geo.material_per_region

    Nmat = length(cs.materials)
    materialNames = Vector{String}()
    for n in range(1,Nmat)
        push!(materialNames,cs.materials[n].name)
    end

    # Calculations of the volume and material array

    nVoxels = [sum(xVoxels),sum(yVoxels),sum(zVoxels)]
    vol = zeros(nVoxels[1],nVoxels[2],nVoxels[3])
    mat = zeros(Int64,nVoxels[1],nVoxels[2],nVoxels[3])
    Δx = zeros(nVoxels[1])
    Δy = zeros(nVoxels[2])
    Δz = zeros(nVoxels[3])

    # Voxels x-width
    i=1
    while i <= Nrx
        if i==1
            rx = range(1,xVoxels[i])
        else
            rx = range(sum(xVoxels[1:i-1])+1,sum(xVoxels[1:i-1])+xVoxels[i])
        end
        for ix in rx
            Δx[ix] = (xB[i+1]-xB[i])/xVoxels[i]
        end
        i += 1
    end

    # Voxels y-width
    j=1
    while j <= Nry
        if j==1
            ry = range(1,yVoxels[j])
        else
            ry = range(sum(yVoxels[1:j-1])+1,sum(yVoxels[1:j-1])+yVoxels[j])
        end
        for iy in ry
            Δy[iy] = (yB[j+1]-yB[j])/yVoxels[j]
        end
        j += 1
    end

    # Voxels z-width
    k=1
    while k <= Nrz
        if k==1
            rz = range(1,zVoxels[k])
        else
            rz = range(sum(zVoxels[1:k-1])+1,sum(zVoxels[1:k-1])+zVoxels[k])
        end
        for iz in rz
            Δz[iz] = (yB[k+1]-yB[k])/zVoxels[k]
        end
        k += 1
    end

    # Volumes
    for ix in range(1,nVoxels[1])
        for iy in range(1,nVoxels[2])
            for iz in range(1,nVoxels[3])
                vol[ix,iy,iz] = Δx[ix] * Δy[iy] * Δz[iz]
            end
        end
    end

    # Material matrix
    i=1
    j=1
    k=1
    while i <= Nrx
        if i==1
            rx = range(1,xVoxels[i])
        else
            rx = range(sum(xVoxels[1:i-1])+1,sum(xVoxels[1:i-1])+xVoxels[i])
        end
        while j <= Nry
            if j==1
                ry = range(1,yVoxels[j])
            else
                ry = range(sum(yVoxels[1:j-1])+1,sum(yVoxels[1:j-1])+yVoxels[j])
            end
            while k <= Nrz
                if k==1
                    rz = range(1,zVoxels[k])
                else
                    rz = range(sum(zVoxels[1:k-1])+1,sum(zVoxels[1:k-1])+zVoxels[k])
                end
                for ix in rx
                    for iy in ry
                        for iz in rz
                            index = findfirst(x->x==idMat[i,j,k],materialNames)
                            if isnothing(index) error("No cross section information for material ",idMat[i,j,k],".") end
                            mat[ix,iy,iz] = index
                        end
                    end
                end
                k += 1
            end
            k = 1
            j += 1
        end
        j = 1
        i += 1
    end

    # Calculations of the x-coordinates of each voxels

    k=1
    xp = zeros(Int64,nVoxels[1])
    x = zeros(nVoxels[1])

    for i in range(1,Nrx)
        for j in range(1,xVoxels[i])
            xp[k]=i
            k=k+1
        end
    end

    for ix in range(1,nVoxels[1])
        stepx = (xB[xp[ix]+1]-xB[xp[ix]])/xVoxels[xp[ix]]
        if xp[ix] == 1
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-1) * stepx
        else
            x[ix] = xB[xp[ix]] + 0.5 * stepx + (ix-sum(xVoxels[1:xp[ix]-1])-1) * stepx
        end
    end

    # Calculations of the y-coordinates of each voxels

    k=1
    yp = zeros(Int64,nVoxels[2])
    y = zeros(nVoxels[2])

    for i in range(1,Nry)
        for j in range(1,yVoxels[i])
            yp[k]=i
            k=k+1
        end
    end

    for iy in range(1,nVoxels[2])
        stepy = (yB[yp[iy]+1]-yB[yp[iy]])/yVoxels[yp[iy]]
        if yp[iy] == 1
            y[iy] = yB[yp[iy]] + 0.5 * stepy + (iy-1) * stepy
        else
            y[iy] = yB[yp[iy]] + 0.5 * stepy + (iy-sum(yVoxels[1:yp[iy]-1])-1) * stepy
        end
    end

    # Calculations of the z-coordinates of each voxels

    k=1
    zp = zeros(Int64,nVoxels[3])
    z = zeros(nVoxels[3])

    for i in range(1,Nrz)
        for j in range(1,zVoxels[i])
            zp[k]=i
            k=k+1
        end
    end

    for iz in range(1,nVoxels[3])
        stepz = (zB[zp[iz]+1]-zB[zp[iz]])/zVoxels[zp[iz]]
        if zp[iz] == 1
            z[iz] = zB[zp[iz]] + 0.5 * stepz + (iz-1) * stepz
        else
            z[iz] = zB[zp[iz]] + 0.5 * stepz + (iz-sum(zVoxels[1:zp[iz]-1])-1) * stepz
        end
    end

    # Save data

    geo.voxels_width["x"] = Δx
    geo.voxels_width["y"] = Δy
    geo.voxels_width["z"] = Δz
    geo.voxels_position["x"] = x
    geo.voxels_position["y"] = y
    geo.voxels_position["z"] = z
    geo.material_per_voxel = mat
    geo.volume_per_voxel = vol

    else
        error("Undefined type of geometry.")
    end

    return geo

end