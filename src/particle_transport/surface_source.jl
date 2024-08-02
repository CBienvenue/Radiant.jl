"""
    surface_source(Q::Array{Union{Array{Float64},Float64}},particle::String,source::Surface_Source,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

Prepare the volume source produced by fixed sources for transport calculations.

# Input Argument(s)
- 'Q': source density.
- 'particle::String': particule type.
- 'source::Volume_Source': volume source information.
- 'cross_sections::Cross_Sections': cross-sections information.
- 'geometry::Geometry': geometry information.
- 'discrete_ordinates::Discrete_Ordinates': discrete_ordinates information.

# Output Argument(s)
- 'Q': source density.
- 'norm': normalization factor.

# Reference(s)
N/A

"""
function surface_source(Q::Array{Union{Array{Float64},Float64}},particle::String,source::Surface_Source,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

    if particle ∉ cross_sections.particles error(string("No cross sections available for ",particle," particle.")) end
    index = findfirst(x -> x == particle,cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    if particle != discrete_ordinates.particle error(string("No methods available for ",particle," particle.")) end
    quadrature_type = discrete_ordinates.quadrature_type
    N = discrete_ordinates.quadrature_order 
    norm = 0.0

    if geometry.type == "cartesian"

    #----
    # Cartesian 1D geometry
    #----
    if geometry.dimension == 1

    # Extract quadrature
    Ω,w = quadrature(N,quadrature_type,1)
    if typeof(Ω) != Vector{Float64} μ = Ω[1] else μ = Ω end
    Nd = length(μ)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    end
    Ω = source.direction

    # Compute the weighting factors
    n_source = 0
    for n in range(1,Nd) 
        if maximum(abs(Ω[1] - μ[n])) < eps()
            n_source = n
            break
        end
    end
    Wn = zeros(Nd)
    if n_source == 0
        for n in range(1,Nd) 
            Wn[n] = (1/abs(μ[n]-Ω[1])^3) 
        end
        Wn_max = maximum(Wn)
        for n in range(1,Nd)
            if Wn[n]/Wn_max < 0.1 Wn[n] = 0.0 end
        end
        Wn = Wn./sum(w .* Wn)
    else
        Wn[n_source] = 1/w[n_source]
    end
    
    # Compute the source
    norm = intensity
    for ig in range(1,Ng)
        if energy_groups[ig] == 1 
            for n in range(1,Nd)
                Q[ig,n,p] += intensity * Wn[n]
            end
        end
    end

    #=
    # Extract quadrature and initialize the surface source array of intensities
    μ,w = quadrature(N,quadrature_type,1)
    number_of_directions = length(μ)
    #Q = zeros(Ng,number_of_directions,2)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    end
    direction = source.direction

    # Find the optimal discrete ordinate
    Δθ = zeros(number_of_directions)
    for n in range(1,number_of_directions)
        Δθ[n] = (μ[n]-direction[1])^2
    end
    if any(x -> abs(x) ≤ 1e-5,Δθ) 
        Np = 1
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        nm[1] = findfirst(x -> abs(x) ≤ eps(),Δθ)
        Pm[1] = 1
    else
        Np = 0
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        P = 1 ./Δθ
        P_max = maximum(P)
        for n in range(1,number_of_directions)
            if P[n]/P_max > 0.05
                Np += 1
                push!(nm,n)
                push!(Pm,P[n]^3)
            end
        end
        Pm = Pm./sum(Pm)
    end

    # Compute the source
    norm = intensity
    for ig in range(1,Ng)
        if energy_groups[ig] == 1 
            for n in range(1,Np)
                Q[ig,nm[n],p] += intensity/w[nm[n]] * Pm[n] 
            end
        end
    end
    =#

    #----
    # Cartesian 2D geometry
    #----
    elseif geometry.dimension == 2

    Ω,w = quadrature(N,quadrature_type,2)
    μ = Ω[1]; η = Ω[2]
    Nd = length(w)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    elseif source.location == "Y-"
        p = 3
    elseif source.location == "Y+"
        p = 4
    end
    Ω₀ = source.direction

    # Compute the weighting factors
    n_source = 0
    for n in range(1,Nd) 
        if maximum([abs(Ω[1][n] - Ω₀[1]),abs(Ω[2][n] - Ω₀[2]),abs(Ω[3][n] - Ω₀[3])]) < eps()
            n_source = n
            break
        end
    end
    Wn = zeros(Nd)
    if n_source == 0
        for n in range(1,Nd) 
            Wn[n] = 1/(abs(Ω[1][n]-Ω₀[1])^3 + abs(Ω[2][n]-Ω₀[2])^3 + abs(Ω[3][n]-Ω₀[3])^3)
        end
        Wn_max = maximum(Wn)
        for n in range(1,Nd)
            if Wn[n]/Wn_max < 0.1 Wn[n] = 0.0 end
        end
        Wn = Wn./sum(w .* Wn)
    else
        Wn[n_source] = 1/w[n_source]
    end

    # Compute the source
    norm = 0
    if source.location == "X-" || source.location == "X+"
        y = geometry.voxels_position["y"]
        Δy = geometry.voxels_width["y"]
        Ny = geometry.number_of_voxels["y"]
        ymin = source.boundaries["y"][1]
        ymax = source.boundaries["y"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Nd)
                if abs(Wn[n]) < eps() continue end
                if Q[ig,n,p] == 0 Q[ig,n,p] = zeros(Ny) end
                xsource = zeros(Ny)
                for iy in range(1,Ny)
                    if ymin <= y[iy] && ymax >= y[iy]
                        xsource[iy] += intensity * Wn[n] * Δy[iy]
                        norm = norm + intensity * Δy[iy]
                    end
                end
                Q[ig,n,p] += xsource
            end
        end
    elseif source.location == "Y-" || source.location == "Y+"
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Nd)
                if abs(Wn[n]) < eps() continue end
                if Q[ig,n,p] == 0 Q[ig,n,p] = zeros(Nx) end
                ysource = zeros(Nx)
                for ix in range(1,Nx)
                    if xmin <= x[ix] && xmax >= x[ix]
                        ysource[ix] += intensity * Wn[n] * Δx[ix]
                        norm = norm + intensity * Δx[ix]
                    end
                end
                Q[ig,nm[n],p] += ysource
            end
        end
    end


    #=
    # Extract quadrature and initialize the surface source array of intensities
    Ω,w = quadrature(N,quadrature_type,2)
    μ = Ω[1]
    η = Ω[2]
    number_of_directions = length(μ)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    elseif source.location == "Y-"
        p = 3
    elseif source.location == "Y+"
        p = 4
    end
    direction = source.direction

    # Finding the optimal discrete ordinates
    Δθ = zeros(number_of_directions)
    for n in range(1,number_of_directions)
        Δθ[n] = (μ[n]-direction[1])^2 + (η[n]-direction[2])^2
    end
    if any(x -> abs(x) ≤ 1e-5,Δθ) 
        Np = 1
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        nm[1] = findfirst(x -> abs(x) ≤ eps(),Δθ)
        Pm[1] = 1
    else
        Np = 0
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        P = 1 ./Δθ
        P_max = maximum(P)
        for n in range(1,number_of_directions)
            if P[n]/P_max > 0.05
                Np += 1
                push!(nm,n)
                push!(Pm,P[n]^3)
            end
        end
        Pm = Pm./sum(Pm)
    end

    # Compute the source
    norm = 0
    if source.location == "X-" || source.location == "X+"
        y = geometry.voxels_position["y"]
        Δy = geometry.voxels_width["y"]
        Ny = geometry.number_of_voxels["y"]
        ymin = source.boundaries["y"][1]
        ymax = source.boundaries["y"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Np)
                if Q[ig,nm[n],p] == 0 Q[ig,nm[n],p] = zeros(Ny) end
                xsource = zeros(Ny)
                for iy in range(1,Ny)
                    if ymin <= y[iy] && ymax >= y[iy]
                        xsource[iy] += intensity/w[nm[n]] * Pm[n] * Δy[iy]
                        norm = norm + intensity * Δy[iy] * Pm[n]
                    end
                end
                Q[ig,nm[n],p] += xsource
            end
        end
    elseif source.location == "Y-" || source.location == "Y+"
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Np)
                if Q[ig,nm[n],p] == 0 Q[ig,nm[n],p] = zeros(Nx) end
                ysource = zeros(Nx)
                for ix in range(1,Nx)
                    if xmin <= x[ix] && xmax >= x[ix]
                        ysource[ix] += intensity/w[nm[n]] * Pm[n] * Δx[ix]
                        norm = norm + intensity * Δx[ix] * Pm[n]
                    end
                end
                Q[ig,nm[n],p] += ysource
            end
        end
    end
    =#

    #----
    # Cartesian 3D geometry
    #----
    elseif geometry.dimension == 3

    Ω,w = quadrature(N,quadrature_type,3)
    μ = Ω[1]; η = Ω[2]; ξ = Ω[3]
    Nd = length(w)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    elseif source.location == "Y-"
        p = 3
    elseif source.location == "Y+"
        p = 4
    elseif source.location == "Z-"
        p = 5
    elseif source.location == "Z+"
        p = 6
    end
    Ω₀ = source.direction

    # Compute the weighting factors
    n_source = 0
    for n in range(1,Nd) 
        if maximum([abs(Ω[1][n] - Ω₀[1]),abs(Ω[2][n] - Ω₀[2]),abs(Ω[3][n] - Ω₀[3])]) < eps()
            n_source = n
            break
        end
    end
    Wn = zeros(Nd)
    if n_source == 0
        for n in range(1,Nd) 
            Wn[n] = 1/(abs(Ω[1][n]-Ω₀[1])^3 + abs(Ω[2][n]-Ω₀[2])^3 + abs(Ω[3][n]-Ω₀[3])^3)
        end
        Wn_max = maximum(Wn)
        for n in range(1,Nd)
            if Wn[n]/Wn_max < 0.1 Wn[n] = 0.0 end
        end
        Wn = Wn./sum(w .* Wn)
    else
        Wn[n_source] = 1/w[n_source]
    end

    # Compute the source
    norm = 0
    if source.location == "X-" || source.location == "X+"
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
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Nd)
                if abs(Wn[n]) < eps() continue end
                if Q[ig,n,p] == 0 Q[ig,n,p] = zeros(Ny,Nz) end
                xsource = zeros(Ny,Nz)
                for iy in range(1,Ny), iz in range(1,Nz)
                    if ymin <= y[iy] && ymax >= y[iy] && zmin <= z[iz] && zmax >= z[iz]
                        xsource[iy,iz] += intensity * Wn[n] * Δy[iy] * Δz[iz]
                        norm = norm + intensity * Δy[iy] * Δz[iz]
                    end
                end
                Q[ig,n,p] += xsource
            end
        end
    elseif source.location == "Y-" || source.location == "Y+"
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        z = geometry.voxels_position["z"]
        Δz = geometry.voxels_width["z"]
        Nz = geometry.number_of_voxels["z"]
        zmin = source.boundaries["z"][1]
        zmax = source.boundaries["z"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Nd)
                if abs(Wn[n]) < eps() continue end
                if Q[ig,n,p] == 0 Q[ig,n,p] = zeros(Nx,Nz) end
                ysource = zeros(Nx,Nz)
                for ix in range(1,Nx), iz in range(1,Nz)
                    if xmin <= x[ix] && xmax >= x[ix] && zmin <= z[iz] && zmax >= z[iz]
                        ysource[ix,iz] += intensity * Wn[n] * Δx[ix] * Δz[iz]
                        norm = norm + intensity * Δx[ix] * Δz[iz]
                    end
                end
                Q[ig,n,p] += ysource
            end
        end
    elseif source.location == "Z-" || source.location == "Z+"
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
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Nd)
                if abs(Wn[n]) < eps() continue end
                if Q[ig,n,p] == 0 Q[ig,n,p] = zeros(Nx,Ny) end
                zsource = zeros(Nx,Ny)
                for ix in range(1,Nx), iy in range(1,Ny)
                    if xmin <= x[ix] && xmax >= x[ix] && ymin <= y[iy] && ymax >= y[iy]
                        zsource[ix,iy] += intensity * Wn[n] * Δx[ix] * Δy[iy]
                        norm = norm + intensity * Δx[ix] * Δy[iy]
                    end
                end
                Q[ig,n,p] += zsource
            end
        end
    end

    #=
    # Extract quadrature and initialize the surface source array of intensities
    Ω,w = quadrature(N,quadrature_type,3)
    μ = Ω[1]
    η = Ω[2]
    ξ = Ω[3]
    number_of_directions = length(μ)

    # Extract source informations
    energy_groups = zeros(Ng)
    for ig in range(1,Ng)
        if ig == source.energy_group energy_groups[ig] = 1 end
    end
    intensity = source.intensity
    if source.location == "X-"
        p = 1
    elseif source.location == "X+"
        p = 2
    elseif source.location == "Y-"
        p = 3
    elseif source.location == "Y+"
        p = 4
    elseif source.location == "Z-"
        p = 5
    elseif source.location == "Z+"
        p = 6
    end
    direction = source.direction

    # Finding the optimal discrete ordinates
    Δθ = zeros(number_of_directions)
    for n in range(1,number_of_directions)
        Δθ[n] = (μ[n]-direction[1])^2 + (η[n]-direction[2])^2 + (ξ[n]-direction[3])^2
    end
    if any(x -> abs(x) ≤ 1e-5,Δθ) 
        Np = 1
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        nm[1] = findfirst(x -> abs(x) ≤ eps(),Δθ)
        Pm[1] = 1
    else
        Np = 0
        nm = zeros(Int64,Np); Pm = zeros(Float64,Np)
        P = 1 ./Δθ
        P_max = maximum(P)
        for n in range(1,number_of_directions)
            if P[n]/P_max > 0.05
                Np += 1
                push!(nm,n)
                push!(Pm,P[n]^3)
            end
        end
        Pm = Pm./sum(Pm)
    end

    # Compute the source
    norm = 0
    if source.location == "X-" || source.location == "X+"
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
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Np)
                if Q[ig,nm[n],p] == 0 Q[ig,nm[n],p] = zeros(Ny,Nz) end
                xsource = zeros(Ny,Nz)
                for iy in range(1,Ny), iz in range(1,Nz)
                    if ymin <= y[iy] && ymax >= y[iy] && zmin <= z[iz] && zmax >= z[iz]
                        xsource[iy,iz] += intensity/w[nm[n]] * Pm[n] * Δy[iy] * Δz[iz]
                        norm = norm + intensity * Δy[iy] * Δz[iz] * Pm[n]
                    end
                end
                Q[ig,nm[n],p] += xsource
            end
        end
    elseif source.location == "Y-" || source.location == "Y+"
        x = geometry.voxels_position["x"]
        Δx = geometry.voxels_width["x"]
        Nx = geometry.number_of_voxels["x"]
        xmin = source.boundaries["x"][1]
        xmax = source.boundaries["x"][2]
        z = geometry.voxels_position["z"]
        Δz = geometry.voxels_width["z"]
        Nz = geometry.number_of_voxels["z"]
        zmin = source.boundaries["z"][1]
        zmax = source.boundaries["z"][2]
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Np)
                if Q[ig,nm[n],p] == 0 Q[ig,nm[n],p] = zeros(Nx,Nz) end
                ysource = zeros(Nx,Nz)
                for ix in range(1,Nx), iz in range(1,Nz)
                    if xmin <= x[ix] && xmax >= x[ix] && zmin <= z[iz] && zmax >= z[iz]
                        ysource[ix,iz] += intensity/w[nm[n]] * Pm[n] * Δx[ix] * Δz[iz]
                        norm = norm + intensity * Δx[ix] * Δz[iz] * Pm[n]
                    end
                end
                Q[ig,nm[n],p] += ysource
            end
        end
    elseif source.location == "Z-" || source.location == "Z+"
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
        for ig in range(1,Ng)
            if energy_groups[ig] != 1 continue end 
            for n in range(1,Np)
                if Q[ig,nm[n],p] == 0 Q[ig,nm[n],p] = zeros(Nx,Ny) end
                zsource = zeros(Nx,Ny)
                for ix in range(1,Nx), iy in range(1,Ny)
                    if xmin <= x[ix] && xmax >= x[ix] && ymin <= y[iy] && ymax >= y[iy]
                        zsource[ix,iy] += intensity/w[nm[n]] * Pm[n] * Δx[ix] * Δy[iy]
                        norm = norm + intensity * Δx[ix] * Δy[iy] * Pm[n]
                    end
                end
                Q[ig,nm[n],p] += zsource
            end
        end
    end

    =#
    else
    error("Cartesian geometries are available only in 1D, 2D and 3D geometries.")
    end
    else
    error("Surface sources are only available in Cartesian geometries.")
    end
    
    return Q,norm

end