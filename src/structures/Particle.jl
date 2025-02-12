abstract type Custom_Particle end
abstract type Photon end
abstract type Electron end
abstract type Positron end

function Photon() Particle(Photon) end
function Electron() Particle(Electron) end
function Positron() Particle(Positron) end

"""
    Particle

Structure used to define a Particle and its properties.

"""
mutable struct Particle

    # Variable(s)
    id      ::Int64
    type    ::Type
    mass    ::Union{Missing,Real} 
    charge  ::Union{Missing,Real}
    spin    ::Union{Missing,Real}

    # Constructor(s)
    function Particle(type::Type)
        this = new()
        this.id = generate_unique_id()
        if type == Photon
            this.mass = 0
            this.charge = 0
            this.spin = 1
        elseif type == Electron
            this.mass = 0.51099895069
            this.charge = -1
            this.spin = 0.5
        elseif type == Positron
            this.mass = 0.51099895069
            this.charge = 1
            this.spin = 0.5
        else
            error("Undefined particle type.")
        end
        this.type = type
        return this
    end
    function Particle(mass::Union{Missing,Real}=missing,charge::Union{Missing,Real}=missing,spin::Union{Missing,Real}=missing)
        this = new()
        this.id = generate_unique_id()
        this.type = Custom_Particle
        this.mass = mass
        this.charge = charge
        this.spin = spin
        return this
    end

end

"""
    get_charge(this::Particle)

Get the charge of the particle.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `charge::Float64` : particle charge.

"""
function get_charge(this::Particle)
    if ismissing(this.charge) error("Particle charge is not defined.") end
    return this.charge
end

"""
    is_photon(this::Particle)

Determines whether the given particle is a photon.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_photon::Bool` : is it a photon or not.

"""
function is_photon(this::Particle)
     if this.type == Photon return true else return false end
end

"""
    is_electron(this::Particle)

Determines whether the given particle is an electron.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_electron::Bool` : is it an electron or not.

"""
function is_electron(this::Particle)
    if this.type == Electron return true else return false end 
end

"""
    is_positron(this::Particle)

Determines whether the given particle is a positron.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_positron::Bool` : is it a positron or not.

"""
function is_positron(this::Particle)
    if this.type == Positron return true else return false end
end

"""
    get_type(this::Particle)

Get the particle type.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `type::Type` : type of particle.

"""
function get_type(this::Particle)
    return this.type
end

"""
    get_id(this::Particle)

Get the particle unique identifier.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `id::Int64` : unique identifier.

"""
function get_id(this::Particle)
    return this.id
end