
mutable struct Material

    # Variable(s)
    name                ::Union{Missing,String}
    density             ::Union{Missing,Float64}
    state_of_matter     ::Union{Missing,String}
    number_of_elements  ::Int64
    elements            ::Vector{String}
    atomic_numbers      ::Vector{Int64}
    weight_fractions    ::Vector{Float64}

    # Constructor(s)
    function Material()
        
        this = new()

        this.name = missing
        this.density = missing
        this.state_of_matter = missing
        this.number_of_elements = 0
        this.elements = Vector{String}()
        this.atomic_numbers = Vector{Int64}()
        this.weight_fractions = Vector{Float64}()

        return this
    end
end

# Method(s)
Base.propertynames(::Material) = 
(
    fieldnames(Material)...,
    :set_name,
    :set_density,
    :set_state_of_matter,
    :add_element,
    :get_density,
    :get_atomic_numbers,
    :weight_fractions
    :get_name
)

function println(this::Material)
    entries = ["Name","Density [g/cm³]","State of matter","Number of elements","Elements","Atomic numbers","Weight fractions"]
    values = [this.name,this.density,this.state_of_matter,this.number_of_elements,this.elements,this.atomic_numbers,this.weight_fractions]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Material")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function set_name(this::Material,name::String)
    this.name = name
end

function set_density(this::Material,density::Float64)
    if density < 0 error("The density should be positive.") end
    this.density = density
end

function set_state_of_matter(this::Material,state_of_matter::String)
    if state_of_matter ∉ ["solid","liquid","gaz"] error("The state of matter is either solid, liquid or gaz.") end
    this.state_of_matter = state_of_matter
end

function add_element(this::Material,symbol::String,weight_fraction::Float64)
    if ~(0 ≤ weight_fraction ≤ 1) error("Weight fraction should have values between 0 and 1.") end
    push!(this.elements,symbol)
    push!(this.atomic_numbers,atomic_number(lowercase(symbol)))
    push!(this.weight_fractions,weight_fraction)
    this.number_of_elements += 1
end

function get_density(this::Material)
    if ismissing(density) error("Unknown density.") end
    return this.density
end

function get_number_of_elements(this::Material)
    return this.number_of_elements
end

function get_atomic_numbers(this::Material)
    return this.atomic_numbers
end

function get_weight_fractions(this::Material)
    return this.weight_fractions
end

function get_name(this::Material)
    return this.name
end

function get_state_of_matter(this::Material)
    return this.state_of_matter
end