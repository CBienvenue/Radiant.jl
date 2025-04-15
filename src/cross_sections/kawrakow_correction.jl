"""
    kawrakow_correction(Z::Int64,Ecutoff::Float64,η::Float64,particle::Particle,
    ai::Vector{Float64}=[],elastic_model::String="rutherford",is_subshells=true)

Gives the Kawrakow correction to the elastic lepton cross-sections.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : kinetic energy of the incoming particle.
- `Ecutoff::Float64` : cutoff energy.
- `η::Float64` : Molière screening parameter.
- `particle::Particle` : particle.
- `ai::Vector{Float64}` : factors for Mott cross-sections.
- `elastic_model::String` : elastic model.
- `is_subshells::Bool` : boolean to enable or not subshell-dependant treatment of the
  inelastic scattering.

# Output Argument(s)
- `ξ::Float64` : Kawrakow correction.

# Reference(s)
- Kawrakow (1997), Improved modeling of multiple scattering in the Voxel Monte Carlo model.
- Karwakow et al. (2021), The EGSnrc Code System: Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function kawrakow_correction(Z::Int64,Ei::Float64,Ecutoff::Float64,η::Float64,particle::Particle,ai::Vector{Float64}=[],elastic_model::String="rutherford",is_subshells::Bool=true)

    #----
    # Initialization
    #----
    Nshells,Zi,Ui,_,_,_ = electron_subshells(Z,~is_subshells)

    #----
    # Compute gM
    #----
    gM = 0
    for δi in range(1,Nshells)
        Wc = Ecutoff # Knock-on production cutoff with Møller or Bhabha
        if is_electron(particle)
            Wmax = (Ei-Ui[δi])/2
            Wmin = Wc
            if Wmin < Wmax
                Jm(x) = ((Ei+2)*((Ei+2)*(x+Ui[δi])*log(x+Ui[δi])-(Ei+2)*(x+Ui[δi])*log(Ei-x+2)+Ui[δi]^2+(Ei+2)*Ui[δi]))/((Ui[δi]+Ei+2)^2*(x+Ui[δi])) + ((Ei+2)^2*log(Ei-x)-(Ei+2)^2*log(Ei-x+2))/4-(Ei*(Ei+2))/(2*(x-Ei)) -((Ei+2)*((Ei+2)*log(Ei-x+2)+x))/(Ei+1)^2 + ((Ei+2)*(2*Ei+1)*(2*Ui[δi]*log(x+Ui[δi])+Ei*(Ui[δi]+Ei+2)*log(Ei-x)-(Ei+2)*(Ui[δi]+Ei)*log(Ei-x+2)))/(2*(Ei+1)^2*(Ui[δi]+Ei)*(Ui[δi]+Ei+2))
                gM += Zi[δi]/Z * (Jm(Wmax) - Jm(Wmin))
            end
        elseif is_positron(particle)
            γ = Ei+1
            b = ((γ-1)/γ)^2
            b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
            b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
            b3 = b * (2*(γ-1)*γ)/(γ+1)^2
            b4 = b * (γ-1)^2/(γ+1)^2
            Wmax = Ei-Ui[δi]
            Wmin = Wc
            if Wmin < Wmax
                Jb(x) = log(x/(Ei-x+2)) + (Ei+2)/Ei * (b1-(Ei+2)/Ei*b2+(Ei+2)^2/Ei^2*b3-(Ei+2)^3/Ei^3*b4)*log(Ei-x+2) -(Ei+2)/Ei^2 * (b2-(Ei+2)/Ei*b3+(Ei+2)^2/Ei^2*b4)*x + (Ei+2)/(2*Ei^3) * (b3-(Ei+2)/Ei*b4)*x^2 -(Ei+2)/(3*Ei^4) * b4*x^3
                gM += Zi[δi]/Z * (Jb(Wmax) - Jb(Wmin))
            end
        else
            error("Unkown particle.")
        end
    end

    #----
    # Compute gR
    #----
    if elastic_model == "boschini"
        gR = 1/3 * ((30*(16*ai[5]*η^3+12*ai[5]*η^2-6*ai[3]*η^2-4*ai[3]*η+2*ai[1]*η+ai[1])*log(2*(η+1))-30*(16*ai[5]*η^3+12*ai[5]*η^2-6*ai[3]*η^2-4*ai[3]*η+2*ai[1]*η+ai[1])*log(2*η)+15*2^(3/2)*atan(1/sqrt(η))*sqrt(η)*(14*ai[4]*η^2+10*ai[4]*η-5*ai[2]*η-3*ai[2])-15*(32*ai[5]+7*2^(5/2)*ai[4])*η^2-5*(24*ai[5]+2^(11/2)*ai[4]-36*ai[3]-15*2^(3/2)*ai[2])*η+20*ai[5]+2^(9/2)*ai[4]+30*ai[3]+5*2^(7/2)*ai[2]-60*ai[1])/10)
    elseif elastic_model == "rutherford"
        gR = (1+2*η)*log(1+1/η)-2
    else
        error("Unknown elastic model.")
    end

    #----
    # Compute Karakow correction
    #----
    ξ = 1 - gM/gR
    return ξ
end