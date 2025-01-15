abstract type Custom_Particle end
abstract type Photon end
abstract type Electron end
abstract type Positron end

function Photon() Particle(Photon) end
function Electron() Particle(Electron) end
function Positron() Particle(Positron) end

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

function get_charge(this::Particle)
    if ismissing(this.charge) error("Particle charge is not defined.") end
    return this.charge
end

function is_photon(this::Particle)
     if this.type == Photon return true else return false end
end
function is_electron(this::Particle)
    if this.type == Electron return true else return false end 
end
function is_positron(this::Particle)
    if this.type == Positron return true else return false end
end
function get_type(this::Particle)
    return this.type
end
function get_id(this::Particle)
    return this.id
end

#=


abstract type Particle end

mutable struct Photon <: Particle
    mass    ::Real 
    charge  ::Real
    spin    ::Real
    function Photon()
        this = new()
        this.mass = 0
        this.charge = 0
        this.spin = 1
        return this
    end
end

mutable struct Electron <: Particle
    mass    ::Real 
    charge  ::Real
    spin    ::Real
    function Electron()
        this = new()
        this.mass = 0.51099895069
        this.charge = -1
        this.spin = 0.5
        return this
    end
end

mutable struct Positron <: Particle
    mass    ::Real 
    charge  ::Real
    spin    ::Real
    function Positron()
        this = new()
        this.mass = 0.51099895069
        this.charge = 1
        this.spin = 0.5
        return this
    end
end

function get_charge(this::Particle)
    return this.charge
end
=#