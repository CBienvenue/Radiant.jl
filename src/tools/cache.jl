# Cache for the Radiant package
const cache_radiant = Ref{Dict{String,Any}}(Dict())
const unique_id_counter = Ref(0)

"""
    generate_unique_id()

Generate a unique identifier through Radiant package.

# Input Argument(s)
N/A

# Output Argument(s)
- `unique_id::Int64` : unique id.

"""
function generate_unique_id()
    unique_id_counter[] += 1
    return unique_id_counter[]
end

"""
    fast_load(file::String)

Load and save data from ".jld2" file in "data/" to cache, if not already in cache, and then 
return it. 

# Input Argument(s)
- `file::String` :  Name of ".jld2" file in "data/".

# Output Argument(s)
- `data::Any` : data contained in the ".jld2" file.

"""
function fast_load(file::String)
    if file[end-4:end] != ".jld2" error("Loading of .jld2 files in data/ folder only.") end
    # Compute and save in cache, if not already in cache
    if ~haskey(cache_radiant[],file[1:end-5])
        cache_radiant[][file[1:end-5]] = load(joinpath(find_package_root(),"data",file))
    end
    # Extract data from cache
    return cache_radiant[][file[1:end-5]]
end


"""
    fast_spline(x::Vector{Float64},f::Vector{Float64},index::Union{Int64,Vector{Int64}},
    str::String)

Compute a spline and save it to cache, if not already in cache, and then return it.

# Input Argument(s)
- `x::Vector{Float64}` : x data to be interpolated.
- `f::Vector{Float64}` : f data to be interpolated.
- `index::Union{Int64,Vector{Int64}}` : index (or list of index) in which the spline is
  kept.
- `str::String` : unique identifier for the splines structure in memory. 

# Output Argument(s)
- `spline::Function` : cubic hermite spline.

"""
function fast_spline(x::Vector{Float64},f::Vector{Float64},index::Union{Int64,Vector{Int64}},str::String)

    # If undefined spline field, generate it.
    id = string("splines_",str)
    if ~haskey(cache_radiant[],id)
        cache_radiant[][id] = Dict()
    end

    # Save the interpolation parameters
    data = cache_radiant[][id]
    N_index = length(index)
    for i in range(1,N_index)
        if ~haskey(data,index[i])
            if i == N_index
                data[index[i]] = cubic_hermite_spline(x,f)
            else
                data[index[i]] = Dict()
            end
        end
        data = data[index[i]]
    end

    # Return the spline
    return data
end