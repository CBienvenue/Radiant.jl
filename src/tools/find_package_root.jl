"""
    find_package_root() 

Find the root of Radiant package.

# Input Argument(s)
N/A

# Output Argument(s)
- `dir::String`: absolute location of the Radiant package root.

"""
function find_package_root()
    current_dir = @__DIR__
    while current_dir != "/"
        if isfile(joinpath(current_dir, "Project.toml"))
            return current_dir
        end
        current_dir = dirname(current_dir)
    end
    error("Package root not found.")
end