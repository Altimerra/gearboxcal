# This file contains the various modules written to perform calculations
# related to the project

module Rolling

export Roller, Material, torque, power

struct Roller
    R # Roller radius
    h0 # initial sample thickness
    d # draft
    w # sample width
    N # RPM
end

struct Material
    k # k factor of material
    n # n factor of material
end
#=
function roll(R,h0,d,w,k,n,N)
    L = sqrt(R*d) # contact length
    e = log(h0/(h0-d)) # true strain
    sigAvg = (k*e^n)/(1+n) # average flow stress
    rf = L*w*sigAvg # Rolling force
    power = 2*pi*rf*L*N / 60 # power
    torque = L*rf/2 # torque
    return torque, power
end
=#
function calc(r::Roller,m::Material)
    L = sqrt(r.R*r.d) # contact length
    e = log(r.h0/(r.h0-r.d)) # true strain
    sigAvg = (m.k*e^m.n)/(1+m.n) # average flow stress
    rf = L*r.w*sigAvg # Rolling force
    torque = L*rf/2 # torque
    power = 2*pi*rf*L*r.N / 60 # power
    return torque, power
end

function torque(r::Roller,m::Material)
    calc(r,m)
    return calc(r,m)[1]
end

function power(r::Roller,m::Material)
    calc(r,m)
    return calc(r,m)[2]
end

end

module Gears
using Classes
using Interpolations
export Gear14c,Gear14f,Gear20f,Gear20s
export plv,wt,tpE,tpT,min_teeth,M,y,beam_safety,cv,lambda,lbd,wtM,er,C,wd,ws,static_safety,Q,K,ww,wear_safety,centre_distance,std_mod
export @solfngen, @vec_gen, @solver

#=
abstract type Gear end

mutable struct Gear14c <: Gear # 14.5 deg composite
    P # power transmitted (Watt)
    dia # diameter (m)
    rpm # Revolutions per minute
    G # gear ratio
    cs # service factor
    T # teeth
    mod # module of gear (mm)
    pinion # whether gear is a pinion or not
    sig # Material strength (Allowable static stress) (N/mm2)
    sige # Flexural endurance limit (N/mm2)
    siges # Surface endurance limit (N/mm2)
    E # Young's modulus (N/mm2)
    met::String # method of manufacture: "mach"(machined) or "cast"
end

mutable struct Gear14f <: Gear # 14.5 deg full depth
    P # power transmitted
    dia # diameter
    rpm # Revolutions per minute
    G # gear ratio
    cs # service factor
    T # teeth
    mod # module of gear
    pinion # whether gear is a pinion or not
    sig # Material strength (stress)
    sige # Flexural endurance limit (N/mm@)
    siges # Surface endurance limit (N/mm2)
    E # Young's modulus (N/mm2)
    met::String # method of manufacture: "mach"(machined) or "cast"
end
=#

@class mutable Gear begin # 20 deg full depth
    P::Real # power transmitted
    dia::Union{Real,Tuple} # diameter
    rpm::Real # Revolutions per minute
    G::Real # gear ratio
    cs::Real # service factor
    T::Real # teeth
    mod::Real # module of gear
    pinion::Bool # whether gear is a pinion or not
    sig::Union{Real,Tuple} # Material strength (stress)
    sige::Union{Real,Tuple} # Flexural endurance limit (N/mm@)
    siges::Union{Real,Tuple} # Surface endurance limit (N/mm2)
    E::Real # Young's modulus (N/mm2)
    met::String # method of manufacture: "mach"(machined) or "cast"
end

@class mutable Gear14c <: Gear begin
end
@class mutable Gear14f <: Gear begin
end
@class mutable Gear20f <: Gear begin
end
@class mutable Gear20s <: Gear begin
end

#=
mutable struct Gear20s <: Gear # 20 deg stub
    P # power transmitted
    dia # diameter
    rpm # Revolutions per minute
    G # gear ratio
    cs # service factor
    T # teeth
    mod # module of gear
    pinion # whether gear is a pinion or not
    sig # Material strength (stress)
    sige # Flexural endurance limit (N/mm@)
    siges # Surface endurance limit (N/mm2)
    E # Young's modulus (N/mm2)
    met::String # method of manufacture: "mach"(machined) or "cast"
end
=#

# pitch line velocity
plv(g::AbstractGear) = (pi*g.dia*g.rpm)/60

# tangential tooth load
wt(g::AbstractGear) = g.P / plv(g) * g.cs

# Minimum pinion teeth (equation)
tpE(g::Gear14c) = 2 / (g.G*(sqrt(1 + (1/g.G + 2)/g.G*sin(deg2rad(14.5))^2)-1))
tpE(g::Gear14f) = 2 / (g.G*(sqrt(1 + (1/g.G + 2)/g.G*sin(deg2rad(14.5))^2)-1))
tpE(g::Gear20f) = 2 / (g.G*(sqrt(1 + (1/g.G + 2)/g.G*sin(deg2rad(20))^2)-1))
tpE(g::Gear20s) = 2 / (g.G*(sqrt(1 + (1/g.G + 2)/g.G*sin(deg2rad(20))^2)-1))

min_teeth(g::AbstractGear) = ceil(tpE(g))

# Minimum pinion teeth (table)
tpT(g::Gear14c) = 12
tpT(g::Gear14f) = 32
tpT(g::Gear20f) = 18
tpT(g::Gear20s) = 14

# Module
M(g::AbstractGear) = g.dia*1000/g.T
# Diameter
diameter(g::AbstractGear) = g.mod*g.T/1000

# TODO insert standard module selector function
std_mod_vals = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.25,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,32,40,50] # ISO series 1
std_mod(g::AbstractGear) = filter(x->x>g.mod,std_mod_vals)[1]

# Lewis form factor
y(g::Gear14c) = 0.124 - 0.684/g.T
y(g::Gear14f) = 0.124 - 0.684/g.T
y(g::Gear20f) = 0.154 - 0.912/g.T
y(g::Gear20s) = 0.175 - 0.841/g.T

# Velocity factor (entire application is low velocity)
cv(g::AbstractGear) = 3/(3 + plv(g))

# b(face width) = lambda*m
# x is between 0,1 and it says where in the range the value is taken
function lambda(g::AbstractGear,x)
if g.met == "mach"
    range = (9.5,12.5)
elseif g.met == "cast"
    range = (6.5,9.5)
else
    range = nothing
end
return (range[2]-range[1])*x+range[1]
end

lbd(g::AbstractGear) = lambda(g,0.5)

# Maximum tooth load (N) (Face width taken as average)
wtM(g::AbstractGear) = cv(g)*g.sig*lbd(g)*(g.mod)^2*pi*y(g)

beam_safety(g::AbstractGear) = wtM(g) .> wt(g)

# Tooth error in action (mm)
err_modules = 4:1:10
err_errors = [0.051,0.055,0.065,0.071,0.078,0.085,0.089]
err_interpol = linear_interpolation(err_modules,err_errors)
er(g::AbstractGear) = err_interpol(g.mod)

# Deformation factor (assuming both gears are made from the same material)
C(g::Gear14c) = 0.107*er(g)*g.E / 2
C(g::Gear14f) = 0.107*er(g)*g.E / 2
C(g::Gear20f) = 0.111*er(g)*g.E / 2
C(g::Gear20s) = 0.115*er(g)*g.E / 2

# Dynamic tooth load (N) (face width taken as average)
wd(g::AbstractGear) = wt(g) + (21*plv(g)*( lbd(g)*g.mod*C(g) + wt(g) ))/(21*plv(g)+sqrt(lbd(g)*g.mod*C(g) + wt(g)))

# static tooth load (face width taken as average)
ws(g::AbstractGear) = g.sige*lbd(g)*(g.mod)^2*pi*y(g)

static_safety(g::AbstractGear) = ws(g) .>= 1.35*wd(g) # Pulsating loads

# Ratio factor (External gears)
Q(g::AbstractGear) = 2*g.G/(g.G + 1)

# Load stress factor (Assuming both gears made of same material)
K(g::Gear14c) = (2*g.siges.^2*sin(deg2rad(14.5)))/(1.4*g.E)
K(g::Gear14f) = (2*g.siges.^2*sin(deg2rad(14.5)))/(1.4*g.E)
K(g::Gear20f) = (2*g.siges.^2*sin(deg2rad(20)))/(1.4*g.E)
K(g::Gear20s) = (2*g.siges.^2*sin(deg2rad(20)))/(1.4*g.E)

# Tooth wear load
ww(g::AbstractGear) = g.dia*1000 * lbd(g)*g.mod*Q(g)*K(g)

wear_safety(g::AbstractGear) = ww(g) .> wd(g)

centre_distance(g::AbstractGear) = g.dia*(1+d.G)/2 # where g is the pinion

macro solfngen(method,sol_for,sol_with)
    return :(
        function val_solve(g,x)
            g.$sol_for = x
            return $method(g) - g.$sol_with
        end
    )
end

macro vec_gen(g,sol_for)
    val = :(g.$sol_for)
    return collect((val-val/10):(val/100):val+val/10)
end

function vec_solve(g::AbstractGear,func,vec::AbstractArray)
    y = [func(g,x) for x in vec]
    min_index = findmin(abs.(y))[2]
    return vec[min_index]
end

end

macro solver(g,func,sol_for) # sol_for is a field that sould contain a range tuple
    ran_tuple = :($g.$sol_for)
    x = :(collect($ran_tuple[1]:(($ran_tuple[2]-$ran_tuple[1])/100):$ran_tuple[2]))
    y = :([$func($g,t) for t in $x])
    min_index = :(findmin(abs.($y))[2])
    return :($x[$min_index])
end
