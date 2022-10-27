# This file contains the various modules written to perform calculations
# related to the project
#=
$$$$$$$\            $$\ $$\ $$\
$$  __$$\           $$ |$$ |\__|
$$ |  $$ | $$$$$$\  $$ |$$ |$$\ $$$$$$$\   $$$$$$\
$$$$$$$  |$$  __$$\ $$ |$$ |$$ |$$  __$$\ $$  __$$\
$$  __$$< $$ /  $$ |$$ |$$ |$$ |$$ |  $$ |$$ /  $$ |
$$ |  $$ |$$ |  $$ |$$ |$$ |$$ |$$ |  $$ |$$ |  $$ |
$$ |  $$ |\$$$$$$  |$$ |$$ |$$ |$$ |  $$ |\$$$$$$$ |
\__|  \__| \______/ \__|\__|\__|\__|  \__| \____$$ |
                                          $$\   $$ |
                                          \$$$$$$  |
                                           \______/
=#
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
#=
 $$$$$$\
$$  __$$\
$$ /  \__| $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$$\
$$ |$$$$\ $$  __$$\  \____$$\ $$  __$$\ $$  _____|
$$ |\_$$ |$$$$$$$$ | $$$$$$$ |$$ |  \__|\$$$$$$\
$$ |  $$ |$$   ____|$$  __$$ |$$ |       \____$$\
\$$$$$$  |\$$$$$$$\ \$$$$$$$ |$$ |      $$$$$$$  |
 \______/  \_______| \_______|\__|      \_______/
=#
module Gears
using Classes
using Interpolations
export Gear14c,Gear14f,Gear20f,Gear20s
export calculate
export centre_distance,gear_ratio

@class mutable Gear begin # 20 deg full depth
    P::Real # power transmitted
    dia::Union{Real,Tuple} # diameter
    rpm::Real # Revolutions per minute
    G::Real # gear ratio
    cs::Real # service factor
    T::Real # teeth
    mod::Real # module of gear
    pinion::Bool # whether gear is a pinion or not
    sig::Union{Real,Tuple,Vector} # Material strength (stress)
    sige::Union{Real,Tuple,Vector} # Flexural endurance limit (N/mm@)
    siges::Union{Real,Tuple,Vector} # Surface endurance limit (N/mm2)
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

lbd(g::AbstractGear) = lambda(g,1)

# Maximum tooth load (N)
wtM(g::AbstractGear) = cv(g)*g.sig*lbd(g)*(g.mod)^2*pi*y(g)

beam_safety(g::AbstractGear) = wtM(g) .> wt(g)

# Tooth error in action (mm)
err_modules = vcat(collect(4:1:10),collect(12:2:20))
err_errors = [0.051,0.055,0.065,0.071,0.078,0.085,0.089,0.097,0.104,0.110,0.114,0.117]
err_interpol = linear_interpolation(err_modules,err_errors)
er(g::AbstractGear) = (g.mod < 4) ? 0.051 : err_interpol(g.mod)

# Deformation factor (assuming both gears are made from the same material)
C(g::Gear14c) = 0.107*er(g)*g.E / 2
C(g::Gear14f) = 0.107*er(g)*g.E / 2
C(g::Gear20f) = 0.111*er(g)*g.E / 2
C(g::Gear20s) = 0.115*er(g)*g.E / 2

# Dynamic tooth load (N)
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

centre_distance(g::AbstractGear) = (g.pinion) ? g.dia*(1+g.G)/2  : g.dia*(1+g.G)/(2*g.G)# where g is the pinion
gear_ratio(g1::AbstractGear,g2::AbstractGear) = (g1.G/g2.G >= 1) ? g1.G/g2.G : g2.G/g1.G
vecgen(ran::Tuple) =  collect(ran[1]:((ran[2]-ran[1])/1000):ran[2])
vecfind(vec::BitVector) = (findfirst(x->x==1,vec)==nothing) ? println("No suitable values in range") : findfirst(x->x==1,vec)

function calculate(g::AbstractGear)
    g.T = min_teeth(g)
    g.mod = M(g)
    g.mod = std_mod(g)
    g.dia = diameter(g)
    g.sig = vecgen(g.sig)
    g.sig = g.sig[vecfind(beam_safety(g))]
    g.sige = vecgen(g.sige)
    g.sige = g.sige[vecfind(static_safety(g))]
    g.siges = vecgen(g.siges)
    g.siges = g.siges[vecfind(wear_safety(g))]
end

function calculate(g1::AbstractGear,g2::AbstractGear)
    if g1.pinion
        p = g1
        g = g2
    elseif g2.pinion
        p = g2
        g = g1
    else
        return
    end
    g.dia = p.dia*p.G
    g.T = p.T*p.G
    g.mod = p.mod
    g.sig = vecgen(g.sig)
    g.sig = g.sig[vecfind(beam_safety(g))]
    g.sige = vecgen(g.sige)
    g.sige = g.sige[vecfind(static_safety(g))]
    g.siges = vecgen(g.siges)
    g.siges = g.siges[vecfind(wear_safety(g))]
end
end

#=
 $$$$$$\  $$\                  $$$$$$\    $$\
$$  __$$\ $$ |                $$  __$$\   $$ |
$$ /  \__|$$$$$$$\   $$$$$$\  $$ /  \__|$$$$$$\    $$$$$$$\
\$$$$$$\  $$  __$$\  \____$$\ $$$$\     \_$$  _|  $$  _____|
 \____$$\ $$ |  $$ | $$$$$$$ |$$  _|      $$ |    \$$$$$$\
$$\   $$ |$$ |  $$ |$$  __$$ |$$ |        $$ |$$\  \____$$\
\$$$$$$  |$$ |  $$ |\$$$$$$$ |$$ |        \$$$$  |$$$$$$$  |
 \______/ \__|  \__| \_______|\__|         \____/ \_______/
=#
module Shafts
using Roots
using QuadGK
using Plots
export Shaft,Force,vecsolve,shear,moment
export calculate

struct Force
    x::Real# X component
    y::Real # Y component
    l::Real # distance from axis origin
end

mutable struct Shaft
    Power::Real # Transmitted power
    N::Real # RPM
    k::Real # di/do diameter ratio
    L::Real # Length
    tau::Real # shaft shear stress
    sig::Real # shaft direct stress
    dia::Union{Tuple,Nothing} # outer diameter
    bear::Tuple # distance of bearings from origin
    forces::Vector{Force} # A list of known forces acting on the shaft
end


Tmax(tau::Real,d_o::Real,k::Real) = pi/16*tau*d_o^3*(1-k^4) # Max torque for a shaft
T(P::Real,N::Real) = P*60/(2*pi*N) # Torque on a shaft
Mmax(sig::Real,d_o::Real,k::Real) = pi/32*sig*d_o^3*(1-k^4) # Max bending moment for a shaft

function vecsolve(s::Shaft)
    b1 = sum([v.x for v in s.forces])
    b2 = sum([v.y for v in s.forces])
    b3 = sum([v.x*v.l for v in s.forces])
    b4 = sum([v.y*v.l for v in s.forces])
    b = [b1;b2;b3;b4]
    r1l = s.bear[1]
    r2l = s.bear[2]
    A = [1 0 1 0;0 1 0 1; r1l 0 r2l 0;0 r1l 0 r2l]
    x = A\(-b)
    r1 = Force(x[1],x[2],r1l)
    r2 = Force(x[3],x[4],r2l)
    return r1,r2
end

function interval(t::Real, a::Real, b::Real)
    heaviside(t) =  0.5 * (sign(t) + 1)
    return heaviside(t-a) - heaviside(t-b)
end

function shear(forces::Vector{Force})
    sort!(forces,by = x -> x.l)
    vx(d::Real,L::Real) = sum([force.x*interval(d,force.l,L) for force in forces])
    vy(d::Real,L::Real) = sum([force.y*interval(d,force.l,L) for force in forces])
    return vx,vy
end

function moment(v::Function,l::Real)
    m(x::Real) = quadgk(z->v(z,l),0,x,rtol=1e-3,atol=1e-3)[1]
    return m
end

Te(m::Real,t::Real,km=2.5,kt=2.25) = sqrt((km*m)^2 + (kt*t)^2) # Guest's theory, related to Tmax
Me(m::Real,t::Real,km=2.5,kt=2.25) = 1/2*(km*m + sqrt((km*m)^2 + (kt*t)^2)) # Rankine's theory, related to Mmax

function calculate(s::Shaft)
    t = T(s.Power,s.N)

    r1,r2 = vecsolve(s)
    append!(s.forces,[r1,r2])
    vx,vy = shear(s.forces)
    mx = moment(vx,s.L)
    my = moment(vy,s.L)
    x = collect(0:(s.L/1000):s.L)
    plvx = plot(x,vx.(x,s.L))
    plvy = plot(x,vy.(x,s.L))
    plmx = plot(x,mx.(x))
    plmy = plot(x,my.(x))
    savefig(plvx,"shear-x.png")
    savefig(plvy,"shear-y.png")
    savefig(plmx,"moment-x.png")
    savefig(plmy,"moment-y.png")
    
    m = findmax( ( findmax(abs.(mx.(x)))[1] , findmax(abs.(my.(x)))[1] ) )[1]
    te = Te(m,t)
    me = Me(m,t)

    t_dia = fzero(x->Tmax(s.tau,x,s.k)-te,0)
    m_dia = fzero(x->Mmax(s.sig,x,s.k)-me,0)
    s.dia = (t_dia,m_dia)
end
end

module Bearings

struct Bearing
    "Radial load"
    wr::Real
    "Axial load"
    wa::Real
    "RPM"
    N::Int
    "Life in hours"
    lh::Int
    "Internal diameter"
    d::Union{Real,Nothing}
    "Static load rating"
    c0::Union{Real,Nothing}
    "Dynamic load rating"
    c::Union{Real,Nothing}
end

"Equivalent radial load"
WeR(wr::Real,wa::Real,x0=0.6,y0=0.5) = x0*wr+y0*wa

"Static load rating for radial ball bearing"
C0(n::Real,Z::Real,D::Real,alpha::Real=0,f_0::Real=12.3,) = f_0*n*Z*D^2*cos(alpha)

#=
"Basic dynamic load for balls smaller than 25.4mm"
C(f_c::Real,n::Real,Z::Real,D::Real,alpha::Real=0) = f_c*(n*cos(alpha))^0.7 * Z*(1/3) * D^1.8 
=#

"Dynamic equivalent radial load"
W(wr::Real,wa::Real,v=1,ks=2,x=1,y=0) = (x*v*wr + y*wa)*ks #ks for medium shock loads

"Dynamic load rating"
C(W::Real,L::Real,k::Real=3) = W*(L/1e6)^(1/k)

"Bearing life multiplier for life other than L10"
a(R::Real) = 6.85*(log(1/R))^(1/1.17)

"Life in revs from life in working hours"
L(lh::Int,N::Real) = 60*N*lh

"Calculates the parameters of the bearing"
function calculate(b::Bearing)
    l10 = L(b.lh,b.N)
    w = W(b.wr,b.wa)
    b.c = C(w,l10)
end

end








#=
macro solfngen(method,sol_for,sol_with)
    return :(
        function val_solve(g,x)
            g.$sol_for = x
            return $method(g) - g.$sol_with
        end
    )
end

#=
macro vec_gen(g,sol_for)
    val = :(g.$sol_for)
    return collect((val-val/10):(val/100):val+val/10)
end

function vec_solve(g::AbstractGear,func,vec::AbstractArray)
    y = [func(g,x) for x in vec]
    min_index = findmin(abs.(y))[2]
    return vec[min_index]
end

=#


# This macro finds a solution for a function from within a range
macro solver(g,func,sol_for) # sol_for is a field that sould contain a range tuple
    #gg = :(esc(g))
    ran_tuple = :($g.$sol_for)
    x = :(collect($ran_tuple[1]:(($ran_tuple[2]-$ran_tuple[1])/100):$ran_tuple[2]))
    y = :([$func($g,t) for t in $x])
    min_index = :(findmin(abs.($y))[2])
    return :($x[$min_index])
end
=#
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
