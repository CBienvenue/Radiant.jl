function voronoi_sphere(Ω)

    Nd = length(Ω[1])
    xn = Vector{Vector{Float64}}(undef,Nd); yn = Vector{Vector{Float64}}(undef,Nd); zn = Vector{Vector{Float64}}(undef,Nd)
    index = Vector{Vector{Vector{Int64}}}(undef,Nd)

    # Function to create a unique key from a float value
    function create_key(x::Float64; precision::Int)
        scale = 10^precision
        return round(Int64, x * scale)
    end

    # Function to create a unique key tuple for the dictionary
    function create_tuple_key(x::Float64, y::Float64, z::Float64; precision::Int)
        return (create_key(x, precision=precision), create_key(y, precision=precision), create_key(z, precision=precision))
    end

    for n in range(1,Nd)

        # Find vertex points
        x_vertex = Vector{Float64}(); y_vertex = Vector{Float64}(); z_vertex = Vector{Float64}()
        pts_index = Vector{Vector{Int64}}()
        j_index = Vector{Int64}(); k_index = Vector{Int64}()
        x1 = Ω[1][n]; y1 = Ω[2][n]; z1 = Ω[3][n]
        for j in range(1,Nd)
            if n == j continue end
            x2 = Ω[1][j]; y2 = Ω[2][j]; z2 = Ω[3][j]
            a1 = x1 - x2; b1 = y1 - y2; c1 = z1 - z2
            d1 = 0.5*(x1^2-x2^2+y1^2-y2^2+z1^2-z2^2)
            for k in range(1,Nd)
                if k == j || k == n continue end
                x3 = Ω[1][k]; y3 = Ω[2][k]; z3 = Ω[3][k]
                a2 = x1 - x3; b2 = y1 - y3; c2 = z1 - z3
                d2 = 0.5*(x1^2-x3^2+y1^2-y3^2+z1^2-z3^2)
                F(x) = [x[1]^2+x[2]^2+x[3]^2-1; a1*x[1]+b1*x[2]+c1*x[3]-d1; a2*x[1]+b2*x[2]+c2*x[3]-d2]
                x⁻ = [x1,y1,z1]
                x = x⁻
                ϵ = Inf
                Nmax = 100
                N = 1
                isinv = true
                while ϵ > 1e-7 && N < Nmax
                    J = [2*x⁻[1] 2*x⁻[2] 2*x⁻[3]; a1 b1 c1; a2 b2 c2]
                    if det(J) != 0
                        inv_J = inv(J)
                    else
                        isinv = false
                        break
                    end
                    x = x⁻ - inv_J*F(x)
                    ϵ = maximum(abs.(F(x)))
                    x⁻ = x
                    N += 1
                end
                if N < Nmax && isinv
                    push!(x_vertex,x[1]); push!(y_vertex,x[2]); push!(z_vertex,x[3])
                    push!(pts_index,[n,j,k])
                else
                    # 3 colinear points, no vertex intersection, skip
                end
            end
        end

        # Remove duplicate
        x_vertex_temp = Vector{Float64}()
        y_vertex_temp = Vector{Float64}()
        z_vertex_temp = Vector{Float64}()
        pts_index_temp = Vector{Vector{Int64}}()
        tuple_dict = Dict{Tuple{Int64, Int64, Int64}, Int64}()
        Nvertex = length(x_vertex)
        precision = 7  # Number of decimal places
        @time for i in 1:Nvertex
            t = create_tuple_key(x_vertex[i], y_vertex[i], z_vertex[i], precision=precision)
            # Check if the tuple key already exists in the dictionary
            if haskey(tuple_dict, t)
                # If exists, update the corresponding pts_index_temp
                index_in_t = tuple_dict[t]
                pts_index_temp[index_in_t] = union(pts_index_temp[index_in_t], pts_index[i])
            else
                # If not exists, add the tuple key to the dictionary and update the temp arrays
                push!(x_vertex_temp, x_vertex[i])
                push!(y_vertex_temp, y_vertex[i])
                push!(z_vertex_temp, z_vertex[i])
                push!(pts_index_temp, pts_index[i])
                tuple_dict[t] = length(x_vertex_temp)
            end
        end
        x_vertex = x_vertex_temp
        y_vertex = y_vertex_temp
        z_vertex = z_vertex_temp
        pts_index = pts_index_temp

        # Find the set of vertex forming the polygon around each quadrature point
        x_vertex_temp = Vector{Float64}(); y_vertex_temp = Vector{Float64}(); z_vertex_temp = Vector{Float64}()
        pts_index_temp = Vector{Vector{Int64}}()
        Nvertex = length(x_vertex)
        for i in range(1,Nvertex)
            d = sqrt.((Ω[1].-x_vertex[i]).^2+(Ω[2].-y_vertex[i]).^2+(Ω[3].-z_vertex[i]).^2)
            if abs(d[n] - minimum(d)) < 1e-7
                push!(x_vertex_temp,x_vertex[i]); push!(y_vertex_temp,y_vertex[i]); push!(z_vertex_temp,z_vertex[i])
                push!(pts_index_temp,sort(pts_index[i]))
            end
        end
        x_vertex = x_vertex_temp; y_vertex = y_vertex_temp; z_vertex = z_vertex_temp
        pts_index = pts_index_temp

        # Order the polygon values
        x_vertex_temp = Vector{Float64}(); y_vertex_temp = Vector{Float64}(); z_vertex_temp = Vector{Float64}()
        pts_index_temp = Vector{Vector{Int64}}()
        Nvertex = length(x_vertex)
        push!(x_vertex_temp,x_vertex[1]); push!(y_vertex_temp,y_vertex[1]); push!(z_vertex_temp,z_vertex[1])
        push!(pts_index_temp,pts_index[1])
        for i in range(2,Nvertex)
            for j in range(1,Nvertex)
                p1 = filter(x -> x != n, pts_index_temp[i-1])
                p2 = filter(x -> x != n, pts_index[j])
                if length(intersect(p1,p2)) == 1 && pts_index[j] ∉ pts_index_temp
                    push!(x_vertex_temp,x_vertex[j]); push!(y_vertex_temp,y_vertex[j]); push!(z_vertex_temp,z_vertex[j])
                    push!(pts_index_temp,pts_index[j])
                    break
                end
            end
        end
        xn[n] = x_vertex_temp
        yn[n] = y_vertex_temp
        zn[n] = z_vertex_temp
        index[n] = pts_index_temp

    end

    # Store information
    voronoi_data = Vector{Dict{String,Any}}(undef,Nd)
    @inbounds for n in range(1,Nd)
        voronoi_data_i = Dict{String,Any}()
        N_triangle = length(xn[n])
        x1 = Ω[1][n]; y1 = Ω[2][n]; z1 = Ω[3][n]
        p1 = [x1,y1,z1]
        voronoi_data_i["Ni_triangle"] = N_triangle
        voronoi_data_i["xi"] = p1
        voronoi_data_i["Si"] = 0.0
        voronoi_data_i["data_triangle_ij"] = Vector{Dict{String,Any}}(undef,N_triangle)
        for m in range(1,N_triangle)
            voronoi_data_ij = Dict{String,Any}()
            if m == 1
                p2 = [xn[n][end],yn[n][end],zn[n][end]]
                index2 = index[n][end]
            else
                p2 = [xn[n][m-1],yn[n][m-1],zn[n][m-1]]
                index2 = index[n][m-1]
            end
            p3 = [xn[n][m],yn[n][m],zn[n][m]]
            index3 = index[n][m]
            voronoi_data_ij["xij⁺"] = p2
            voronoi_data_ij["xij⁻"] = p3
            voronoi_data_ij["ℓij"] = acos(min(max(dot(p2,p3),-1),1))

            # Arc on unit sphere, spherical law of cosines and Girard's theorem
            a = acos(min(max(dot(p1,p2),-1),1)); b = acos(min(max(dot(p2,p3),-1),1)); c = acos(min(max(dot(p3,p1),-1),1))
            A = acos(min(max((cos(a)-cos(c)*cos(b))/(sin(c)*sin(b)),-1),1))
            B = acos(min(max((cos(b)-cos(a)*cos(c))/(sin(a)*sin(c)),-1),1))
            C = acos(min(max((cos(c)-cos(a)*cos(b))/(sin(a)*sin(b)),-1),1))
            voronoi_data_ij["Sij"] = A + B + C - π
            voronoi_data_i["Si"] += A + B + C - π

            # Find the corresponding node on the other side of the boundary
            is_p4 = false
            index4 = 0
            p4 = Vector{Float64}(undef,3)
            for i in range(1,Nd)
                if i == n continue end
                if index2 in index[i] && index3 in index[i]
                    p4 = [Ω[1][i],Ω[2][i],Ω[3][i]]
                    index4 = i
                    is_p4 = true
                    break
                end
            end
            if ~is_p4 error() end
            voronoi_data_ij["index"] = index4
            voronoi_data_ij["xij"] = p4
            voronoi_data_ij["Δxij"] = acos(min(max(dot(p1,p4),-1),1))
            voronoi_data_i["data_triangle_ij"][m] = voronoi_data_ij
        end
        voronoi_data[n] = voronoi_data_i
    end

    return voronoi_data
end

function fokker_planck_weights_3D(Ω,w,voronoi_data)
    
    Nd = length(w)

    # Mapping for each edge such as (min(n,m),max(n,m)) => k-edge
    edge_list = Vector{Tuple{Int64,Int64}}()
    for n in range(1,Nd)
        Nn_triangle = voronoi_data[n]["Ni_triangle"]
        for j in range(1,Nn_triangle)
            m = voronoi_data[n]["data_triangle_ij"][j]["index"]
            if (min(n,m),max(n,m)) ∉ edge_list
                push!(edge_list,(min(n,m),max(n,m)))
            end
        end
    end

    # Compute matrix system
    Nunk = length(edge_list) 
    Neqs = 3*Nd
    Γ = zeros(Neqs,Nunk)
    Q = zeros(Neqs)
    index_eq = 1
    for n in range(1,Nd), k in range(1,3)

        if k == 0
            R = ones(length(Ω[1]))
            ℓ = 0
        elseif k == 1
            R = Ω[1]
            ℓ = 1
        elseif k == 2
            R = Ω[2]
            ℓ = 1
        elseif k == 3
            R = Ω[3]
            ℓ = 1
        else
            error()
        end

        Nn_triangle = voronoi_data[n]["Ni_triangle"]
        for j in range(1,Nn_triangle)
            m = voronoi_data[n]["data_triangle_ij"][j]["index"]
            index_unk = findfirst(x -> x == (min(n,m),max(n,m)), edge_list)
            Γ[index_eq,index_unk] -= R[n]
            Γ[index_eq,index_unk] += R[m]
            Q[index_eq] = -2*w[n]*R[n]
        end
        index_eq += 1
    end
    γ = pinv(Γ)*Q

    return γ, edge_list

end

function fokker_planck_weights_2D(Ω,w,Ω_3D,w_3D,voronoi_data)

    Nd = length(w)
    Nd_3D = length(w_3D)

    # 2D mapping
    map_3D_to_2D = Vector{Int64}(undef,Nd_3D)
    map_2D_to_3D = Vector{Int64}(undef,Nd)
    for n_3D in range(1,Nd_3D)
        ϵ_ξ = ones(Nd)*Inf
        for j in range(1,Nd)
            ϵ_ξ[j] = abs(Ω_3D[1][n_3D] - Ω[1][j])^2 + abs(Ω_3D[2][n_3D] - Ω[2][j])^2 + abs(Ω_3D[3][n_3D] - sign(Ω_3D[3][n_3D]) * Ω[3][j])^2
        end
        map_3D_to_2D[n_3D] = argmin(ϵ_ξ)
        if Ω_3D[3][n_3D] ≥ 0
            map_2D_to_3D[argmin(ϵ_ξ)] = n_3D
        end
    end

    # Mapping for each edge such as (min(n,m),max(n,m)) => k-edge
    edge_list = Vector{Tuple{Int64,Int64}}()
    for n in range(1,Nd_3D)
        Nn_triangle = voronoi_data[n]["Ni_triangle"]
        if Ω_3D[3][n] < 0 continue end
        for j in range(1,Nn_triangle)
            m = voronoi_data[n]["data_triangle_ij"][j]["index"]
            if Ω_3D[3][m] < 0 continue end
            if (min(map_3D_to_2D[n],map_3D_to_2D[m]),max(map_3D_to_2D[n],map_3D_to_2D[m])) ∉ edge_list
                push!(edge_list,(min(map_3D_to_2D[n],map_3D_to_2D[m]),max(map_3D_to_2D[n],map_3D_to_2D[m])))
            end
        end
    end

    # Compute matrix system
    Nunk = length(edge_list) 
    Neqs = 3*Nd
    Γ = zeros(Neqs,Nunk)
    Q = zeros(Neqs)
    index_eq = 1
    for n in range(1,Nd), k in range(0,2)

        if k == 0
            R = ones(length(Ω[1]))
            ℓ = 0
        elseif k == 1
            R = Ω[1]
            ℓ = 1
        elseif k == 2
            R = Ω[2]
            ℓ = 1
        elseif k == 3
            R = (3*Ω[1].^2 .- 1)./2
            ℓ = 2
        elseif k == 4
            R = sqrt(3) .* Ω[1] .* Ω[2]
            ℓ = 2
        elseif k == 5
            R = sqrt(3)/2 .* (Ω[2].^2 .- Ω[3].^2)
            ℓ = 2
        elseif k == 6
            R = Ω[1] .* (5*Ω[1].^2 .- 3)./2
            ℓ = 3
        elseif k == 7
            R = sqrt(3/8) .* Ω[2] .* (5*Ω[1].^2 .- 1)
            ℓ = 3
        elseif k == 8
            R = sqrt(15)/2 .* Ω[1] .* (Ω[2].^2 .- Ω[3].^2)
            ℓ = 3
        elseif k == 9
            R = sqrt(5/8) .* Ω[2] .* (Ω[2].^2 .- 3*Ω[3].^2)
            ℓ = 3
        elseif k == 10
            R = 35/8 .* Ω[1].^4 - 15/4 .* Ω[1].^2 .+ 3/8
            ℓ = 4
        elseif k == 11
            R = sqrt(10)/4 .* Ω[1] .* Ω[2] .* (7 .* Ω[1].^2 .- 3)
            ℓ = 4
        elseif k == 12
            R = sqrt(5)/4 .* (Ω[1].^2 + 2 .* Ω[2].^2 .- 1) .* (7 .* Ω[1].^2 .- 1)
            ℓ = 4
        elseif k == 13
            R = sqrt(70)/4 .* Ω[1] .* Ω[2] .* (3 .* Ω[1].^2 + 4 .* Ω[2].^2 .- 3)
            ℓ = 4
        elseif k == 14
            R = sqrt(35)/8 .* (Ω[1].^4 + (8 .* Ω[2].^2 .- 2) .* Ω[1].^2 + 8 .* Ω[2].^4 - 8 .* Ω[2].^2 .+ 1)
            ℓ = 4
        else
            error()
        end

        Nn_triangle = voronoi_data[map_2D_to_3D[n]]["Ni_triangle"]
        for j in range(1,Nn_triangle)
            m_3D = voronoi_data[map_2D_to_3D[n]]["data_triangle_ij"][j]["index"]
            if Ω_3D[3][m_3D] < 0 continue end
            m = map_3D_to_2D[m_3D]
            if n == m continue end
            index_unk = findfirst(x -> x == (min(n,m),max(n,m)), edge_list)
            Γ[index_eq,index_unk] -= R[n]
            Γ[index_eq,index_unk] += R[m]
        end

        Q[index_eq] = -ℓ*(ℓ+1)/2 * 2 * w[n] * R[n]
    
        index_eq += 1
    end
    γ = pinv(Γ)*Q

    return γ, edge_list

end