module Radiant

    # External packages
    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using JLD2

    # Find root of Radiant
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

    # Global unique ID generator
    const unique_id_counter = Ref(0)
    function generate_unique_id()
        unique_id_counter[] += 1
        return unique_id_counter[]
    end

    # Function and structures
    include("./structures/Structures.jl")
    include("./tools/Tools.jl")
    include("./cross_sections/Cross_sections.jl")
    include("./particle_transport/Particle_transport.jl")

end
