abstract type Custom_Particle end
abstract type Photon end
abstract type Electron end
abstract type Positron end
abstract type Proton end
abstract type Antiproton end
abstract type Alpha end
abstract type Muon end
abstract type Antimuon end

function Photon() Particle(Photon) end
function Electron() Particle(Electron) end
function Positron() Particle(Positron) end
function Proton() Particle(Proton) end
function Antiproton() Particle(Antiproton) end
function Alpha() Particle(Alpha) end
function Muon() Particle(Muon) end
function Antimuon() Particle(Antimuon) end

function Photon(tag::String) p = Particle(Photon); p.tag = tag; return p end
function Electron(tag::String) p = Particle(Electron); p.tag = tag; return p end
function Positron(tag::String) p = Particle(Positron); p.tag = tag; return p end
function Proton(tag::String) p = Particle(Proton); p.tag = tag; return p end
function Antiproton(tag::String) p = Particle(Antiproton); p.tag = tag; return p end
function Alpha(tag::String) p = Particle(Alpha); p.tag = tag; return p end
function Muon(tag::String) p = Particle(Muon); p.tag = tag; return p end
function Antimuon(tag::String) p = Particle(Antimuon); p.tag = tag; return p end

"""
    Particle

Structure used to define a Particle and its properties.

"""
mutable struct Particle

    # Variable(s)
    tag     ::String
    type    ::Type
    mass    ::Union{Missing,Real}
    charge  ::Union{Missing,Real}

    # Constructor(s)
    function Particle(type::Type)
        this = new()
        if type == Photon
            this.tag = "photon"
            this.mass = 0
            this.charge = 0
        elseif type == Electron
            this.tag = "electron"
            this.mass = 0.51099895069
            this.charge = -1
        elseif type == Positron
            this.tag = "positron"
            this.mass = 0.51099895069
            this.charge = 1
        elseif type == Proton
            this.tag = "proton"
            this.mass = 938.2720894
            this.charge = 1
        elseif type == Antiproton
            this.tag = "antiproton"
            this.mass = 938.2720894
            this.charge = -1
        elseif type == Alpha
            this.tag = "alpha"
            this.mass = 3727.3794118
            this.charge = 2
        elseif type == Muon
            this.tag = "muon"
            this.mass = 105.6583755
            this.charge = -1
        elseif type == Antimuon
            this.tag = "antimuon"
            this.mass = 105.6583755
            this.charge = 1
        else
            error("Undefined particle type.")
        end
        this.type = type
        return this
    end
    function Particle(tag::String,mass::Union{Missing,Real}=missing,charge::Union{Missing,Real}=missing)
        if isempty(tag) error("Particle tag cannot be empty.") end
        this = new()
        this.tag = tag
        this.type = Custom_Particle
        this.mass = mass
        this.charge = charge
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
    get_mass(this::Particle)

Get the mass of the particle [in MeV].

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `charge::Float64` : particle mass.

"""
function get_mass(this::Particle)
    if ismissing(this.mass) error("Particle mass is not defined.") end
    return this.mass
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
    is_proton(this::Particle)

Determines whether the given particle is a proton.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_proton::Bool` : is it a proton or not.

"""
function is_proton(this::Particle)
    if this.type == Proton return true else return false end
end

"""
    is_antiproton(this::Particle)

Determines whether the given particle is an antiproton.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_antiproton::Bool` : is it an antiproton or not.

"""
function is_antiproton(this::Particle)
    if this.type == Antiproton return true else return false end
end

"""
    is_alpha(this::Particle)

Determines whether the given particle is an alpha particle.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_alpha::Bool` : is it an alpha particle or not.

"""
function is_alpha(this::Particle)
    if this.type == Alpha return true else return false end
end

"""
    is_muon(this::Particle)

Determines whether the given particle is a muon.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_muon::Bool` : is it a muon or not.

"""
function is_muon(this::Particle)
    if this.type == Muon return true else return false end
end

"""
    is_antimuon(this::Particle)

Determines whether the given particle is an antimuon.

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `is_antimuon::Bool` : is it a antimuon or not.

"""
function is_antimuon(this::Particle)
    if this.type == Antimuon return true else return false end
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
    get_tag(this::Particle)

Get the particle tag (string identifier).

# Input Argument(s)
- `this::Particle` : particle.

# Output Argument(s)
- `tag::String` : tag.

"""
function get_tag(this::Particle)
    return this.tag
end

"""
    set_tag(this::Particle,tag::String)

Set the particle tag (string identifier).

# Input Argument(s)
- `this::Particle` : particle.
- `tag::String` : tag.

"""
function set_tag(this::Particle,tag::String)
    if isempty(tag) error("Particle tag cannot be empty.") end
    this.tag = tag
end