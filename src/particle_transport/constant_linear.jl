"""
    constant_linear(x::Float64,y::Float64)

Compute weighting factors for linear flux expansion on 2D finite-element.

# Input Argument(s)
- 'x::Float64': centroid x-position of the linear part.
- 'y::Float64': centroid y-position of the linear part.

# Output Argument(s)
- 'Q::Vector{Float64}': vector containing Q-values.
- 'T::Vector{Float64}': vector containing T-values.

# Reference(s)
- Bienvenue (2023) : Adaptive Gradient-Driven Coupled Linear Schemes and their Usefulness
  for Charged Particle Transport.

"""
function constant_linear(x::Float64,y::Float64)

subdomain = which_subdomain(x,y)
δ1,δ2 = compute_δ(x,y,subdomain)
Q,T = compute_QT(δ1,δ2,subdomain)

return Q,T
end

#----
# To compute the T factor corresponding to a (δ1,δ2) pair
#----
function compute_QT(δ1,δ2,subdomain)

    Q = zeros(2); T = zeros(2)

    if subdomain == "A⁺⁺"
        Q[1] += QA⁺(δ1,δ2)
        Q[2] += QA⁺(δ2,δ1)
        T[1] += TA⁺(δ1,δ2)
        T[2] += TA⁺(δ2,δ1)
    elseif subdomain == "A⁺⁻"
        Q[1] += QA⁺(δ1,δ2)
        Q[2] += QA⁻(δ2,δ1)
        T[1] -= TA⁺(δ1,δ2)
        T[2] -= TA⁻(δ2,δ1)
    elseif subdomain == "A⁻⁺"
        Q[1] += QA⁻(δ1,δ2)
        Q[2] += QA⁺(δ2,δ1)
        T[1] -= TA⁻(δ1,δ2)
        T[2] -= TA⁺(δ2,δ1)
    elseif subdomain == "A⁻⁻"
        Q[1] += QA⁻(δ1,δ2)
        Q[2] += QA⁻(δ2,δ1)
        T[1] += TA⁻(δ1,δ2)
        T[2] += TA⁻(δ2,δ1)
    elseif subdomain == "Bx⁺"
        Q[1] += QB⁺(δ1,δ2)
        Q[2] += QB(δ1,δ2)
        T[1] += TB⁺(δ1,δ2)
        T[2] += TB(δ1,δ2) 
    elseif subdomain == "Bx⁻"
        Q[1] += QB⁻(δ1,δ2)
        Q[2] += QB(δ1,δ2)
        T[1] += TB⁻(δ1,δ2)
        T[2] -= TB(δ1,δ2)
    elseif subdomain == "By⁺"
        Q[1] += QB(δ1,δ2)
        Q[2] += QB⁺(δ1,δ2)
        T[1] += TB(δ1,δ2)
        T[2] += TB⁺(δ1,δ2) 
    elseif subdomain == "By⁻"
        Q[1] += QB(δ1,δ2)
        Q[2] += QB⁻(δ1,δ2)
        T[1] -= TB(δ1,δ2)
        T[2] += TB⁻(δ1,δ2)
    elseif subdomain == "C⁺⁺"
        Q[1] += QC⁺(δ1)
        Q[2] += QC⁺(δ2)
        T[1] += TC⁺(δ1,δ2)
        T[2] += TC⁺(δ2,δ1)
    elseif subdomain == "C⁺⁻"
        Q[1] += QC⁺(δ1)
        Q[2] += QC⁻(δ2)
        T[1] -= TC⁺(δ1,δ2)
        T[2] -= TC⁻(δ2,δ1)
    elseif subdomain == "C⁻⁺"
        Q[1] += QC⁻(δ1)
        Q[2] += QC⁺(δ2)
        T[1] -= TC⁻(δ1,δ2)
        T[2] -= TC⁺(δ2,δ1)
    elseif subdomain == "C⁻⁻"
        Q[1] += QC⁻(δ1)
        Q[2] += QC⁻(δ2)
        T[1] += TC⁻(δ1,δ2)
        T[2] += TC⁻(δ2,δ1)
    else
        error("Unknown subdomain: ",subdomain)
    end

    if any(x->isnan(x),Q) || any(x->isinf(x),Q) ||any(x->isnan(x),T)|| any(x->isinf(x),T) || isinf(δ1) || isinf(δ2) || isnan(δ1) || isnan(δ2) error("Values for subdomain ",subdomain," for (Qx,Qy)=(",Q[1],",",Q[2],") and (Tx,Ty)=(",T[1],",",T[2],") and (δ1,δ2)=(",δ1,",",δ2,").") end

    return Q,T

end

QA⁺(δi,δj) = (6-2*δi^2*δj)/(3*δi^3*δj+6*(1-δi^2*δj))
QA⁻(δi,δj) = (2*δi^2*δj+6*(1-δi*δj))/(3*δi^3*δj+6*(1-δi^2*δj))
QB⁺(δi,δj) = if (abs(δi-δj) < 1e-5) return (1+δi)/((1-δi)*(1+2*δi)) else return (1-(δi^2+δj^2+δi*δj)/3)/(1-(δi^2+δj^2+δi*δj)+0.5*(δi^3+δj^3+δi^2*δj+δi*δj^2))       end
QB(δi,δj)  = if (abs(δi-δj) < 1e-5) return 1.0                      else return (1-δi/3-2*δj/3)/(1-0.5*(δi+δj))                                                     end
QB⁻(δi,δj) = if (abs(δi-δj) < 1e-5) return 1/(1+2*δi)               else return (1-δi-δj+(δi^2+δj^2+δi*δj)/3)/(1-(δi^2+δj^2+δi*δj)+0.5*(δi^3+δj^3+δi^2*δj+δi*δj^2)) end
QC⁺(δi)    = 2*(2+δi)/(3*(1-δi^2))
QC⁻(δi)    = 2/(3*(1+δi))

TA⁺(δi,δj) = -δi^2*δj*(δi^3*δj^2+12*δi*δj-30*δi-20*δj+40)/(5*(δi^3*δj-2*δi^2*δj+2)^2)
TA⁻(δi,δj) = -δi*δj*(δi^4*δj^2-8*δi^3*δj^2+(10*δj^2-12*δj+30)*δi^2+(40*δj-80)*δi-40*δj+60)/(5*(δi^3*δj-2*δi^2*δj+2)^2)
TB⁺(δi,δj) = if (abs(δi+δj-2) < 1e-5) return 0.0 else return -(δi-δj)*((δi^4+4*δi^3*δj+10*δi^2*δj^2+4*δi*δj^3+δj^4)-6*(3*δi^2+4*δi*δj+3*δj^2)+20*(δi+δj))/(5*(δi^3+(δj-2)*δi^2+(δj^2-2*δj)*δi+δj^3-2*δj^2+2)^2) end
TB⁻(δi,δj) = if (abs(δi-δj) < 1e-5) return 0.0 else return (δi-δj)*((δi^4+4*δi^3*δj+10*δi^2*δj^2+4*δi*δj^3+δj^4)-8*(δi^3+4*δi^2*δj+4*δi*δj^2+δj^3)+4*(7*δi^2+16*δi*δj+7*δj^2)-40*(δi+δj)+20)/(5*(δi^3+(δj-2)*δi^2+(δj^2-2*δj)*δi+δj^3-2*δj^2+2)^2) end
TB(δi,δj)  = if (abs(δi-δj) < 1e-5) return 0.0 else return -(δi-δj)*(δi^2+6*δi*δj+3*δj^2-8*δi-12*δj+10)/(5*(δi+δj-2)^2) end
TC⁺(δi,δj) = if (abs(δi-δj) < 1e-5) return -(δi^2+6*δi+3)/(5*(1+δi)^2) else return -(δi^2+6*δi+3)*(δj-1)/(5*(1+δi)^2*(δi-1)) end
TC⁻(δi,δj) = -(δi-1)*(δj-1)/(5*(1+δi)^2)

#----
# To compute δ lenghts
#----
function compute_δ(x,y,subdomain)

    # Initialization
    δ1 = 0; δ2 = 0
    
    # Compute δ lengths
    if subdomain == "A⁺⁺"
        δ1,δ2 = δA(x,y)
    elseif subdomain == "A⁺⁻"
        δ1,δ2 = δA(x,-y)
    elseif subdomain == "A⁻⁺"
        δ1,δ2 = δA(-x,y)
    elseif subdomain == "A⁻⁻"
        δ1,δ2 = δA(-x,-y)
    elseif subdomain == "Bx⁺"
        δ1,δ2 = δB(x,y)
    elseif subdomain == "Bx⁻"
        δ1,δ2 = δB(-x,y)
    elseif subdomain == "By⁺"
        δ1,δ2 = δB(y,x)
    elseif subdomain == "By⁻"
        δ1,δ2 = δB(-y,x)
    elseif subdomain == "C⁺⁺"
        δ1,δ2 = δC(x,y)
    elseif subdomain == "C⁺⁻"
        δ1,δ2 = δC(x,-y)
    elseif subdomain == "C⁻⁺"
        δ1,δ2 = δC(-x,y)
    elseif subdomain == "C⁻⁻"
        δ1,δ2 = δC(-x,-y)
    else
        error("Unknown subdomain.")
    end

    if (δ1 > 1 && δ1 < 1.0000001) δ1 = 1 end
    if (δ1 < 0 && δ1 > -0.0000001) δ1 = 0 end
    if (δ2 > 1 && δ2 < 1.0000001) δ2 = 1 end
    if (δ2 < 0 && δ2 > -0.0000001) δ2 = 0 end
    if isnan(δ1) || isnan(δ2) || isinf(δ1) || isinf(δ2) || δ1 < 0 || δ1 > 1 || δ2 < 0 || δ2 > 1 error("Values for subdomain ",subdomain," for (δ1,δ2)=(",δ1,",",δ2,") and (x,y)=(",x,",",y,").") end

    return δ1,δ2

end

function δA(x,y)
    if abs(x) < 1e-10 || abs(y) < 1e-10
        δ1 = 0
        δ2 = 0
    else
        Δ = -48*(9*x^4*y^2+14*x^3*y^3+9*x^2*y^4) - 144*(3*x^4*y+5*x^3*y^2+5*x^2*y^3+3*x*y^4) - 36*(3*x^4-10*x^3*y-18*x^2*y^2-10*x*y^3+3*y^4) + 12*(7*x^3-15*x^2*y-15*x*y^2+7*y^3) - 27*(x^2-2*x*y+y^2)
        τ = 64*x^3*y^3 + 48*(x^3*y^2+x^2*y^3) - 24*(x^3*y+4*x^2*y^2+x*y^3) - 4*(2*x^3-3*x^2*y-3*x*y^2+2*y^3) + 4*x*y*sqrt(Complex(Δ))
        ζ = real((1+im*sqrt(3))*τ^(1/3) - 8*x*y + 4*(4*x^2*y^2+2*x^2*y+2*x*y^2+x^2+y^2-x*y)*(1-im*sqrt(3))*τ^(-1/3))
        δ1 = 1 - (4*x + ζ)/(8*y)
        δ2 = 1 - (4*y + ζ)/(8*x)
    end
    return δ1,δ2
end

function δB(x,y)
    δ1 = (12*y^2+(6-12*x)*y+2*x)/(12*y^2+1)
    δ2 = (12*y^2-(6-12*x)*y+2*x)/(12*y^2+1)
    return δ1,δ2
end

function δC(x,y)
    δ1 = 3*x-1/2
    δ2 = 3*y-1/2
    return δ1,δ2
end

#----
# To determine the subdomain corresponding to a centroid position (x,y) 
#----
function which_subdomain(x,y)

    # Initialization
    if x < -1/2 || x > 1/2 || y < -1/2 || y > 1/2 error("The centroid position of the finite element constant part is out of the finite element domain.") end
    subdomain = ""

    # Find subdomain
    if -1/2 <= x <= -1/6
        if -1/2 <= y <= -1/6
            subdomain = "C⁻⁻"
        elseif -1/6 < y < 1/6
            subdomain = "Bx⁻"
        elseif 1/6 <= y <= 1/2
            subdomain = "C⁻⁺"
        end
    elseif -1/6 < x <= 0
        if -1/2 <= y <= y_out(-x)
            subdomain = "By⁻"
        elseif y_out(-x) < y < -y_in(-x)
            subdomain = "A⁻⁻"
        elseif -y_in(-x) <= y <= y_in(-x)
            subdomain = "Bx⁻"
        elseif y_in(-x) < y < -y_out(-x)
            subdomain = "A⁻⁺"
        elseif -y_out(-x) <= y <= 1/2
            subdomain = "By⁺"
        end
    elseif 0 < x < 1/6
        if -1/2 <= y <= y_out(x)
            subdomain = "By⁻"
        elseif y_out(x) < y < -y_in(x)
            subdomain = "A⁺⁻"
        elseif -y_in(x) <= y <= y_in(x)
            subdomain = "Bx⁺"
        elseif y_in(x) < y < -y_out(x)
            subdomain = "A⁺⁺"
        elseif -y_out(x) <= y <= 1/2
            subdomain = "By⁺"
        end
    elseif 1/6 <= x <= 1/2
        if -1/2 <= y <= -1/6
            subdomain = "C⁺⁻"
        elseif -1/6 < y < 1/6
            subdomain = "Bx⁺"
        elseif 1/6 <= y <= 1/2
            subdomain = "C⁺⁺"
        end
    end

    return subdomain
end

# Boundaries of A-type subdomains
y_in(x) = (-6 * x - 3 + sqrt(36 * x ^ 2 - 60 * x + 9)) / (36 * x - 30 - 6 * sqrt(36 * x ^ 2 - 60 * x + 9))
y_out(x) = (6 * x ^ 2 - 3 * x) / (6 * x + 1)