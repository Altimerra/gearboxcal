begin
import Pkg
Pkg.activate(".")
include("calculator.jl")
import .Rolling, .Gears, .Shafts, .Bearings
end

begin
sheet = Rolling.Roller(30e-3, 6e-3, 1e-3, 50e-3,1*60/5 * 0.15)
wire = Rolling.Roller(30e-3, 8e-3, 1e-3, 8e-3,1*60/5)
brass = Rolling.Material(895e6, 0.49)

torque_sheet = Rolling.torque(sheet,brass)
power_sheet = Rolling.power(sheet,brass)

torque_sheet = Rolling.torque(wire,brass)
power_sheet = Rolling.power(wire,brass)
end

begin
power_motor = 90
cs = 1.25 # Medium shock
pres_angle = deg2rad(20)

gear(
    d::Real,
    N::Real,
    G::Real,
    Pin::Bool
)  =
Gears.Gear20s(
    power_motor, # Power
    d, # diameter
    N, # RPM
    G, # Gear ratio
    cs, # Service factor
    10, # Teeth (nominal)
    1, # module(nominal)
    Pin, # Pinion
    (50,500), # static strength
    (80,700), # Flexural endurance limit
    (600,1050), # Surface endurace
    200000, # Elastic modulus
    "mach", # Manufacturing
    nothing, # face width
    0 # Min face width factor
)


ga = gear(24e-3,60,5,true)

Gears.calculate(ga)
gb = gear(50e-3,60,5,false)

Gears.calculate(ga,gb)

Gears.centre_distance(gb)
end

begin
shaft_len = 50e-3
# Define the forces acting on shaft
f1 = Shafts.Force(Gears.wt(ga),Gears.wt(ga)*tan(pres_angle),shaft_len/3)

s1 = Shafts.Shaft(power_motor,ga.rpm,0.7,shaft_len,50e6,100e6,nothing,(0,shaft_len),[f1])
Shafts.calculate(s1)
end

s1
#=
begin
f1 = Shafts.Force(2,0,1)
f2 = Shafts.Force(-2,0,2)
s1 = Shafts.Shaft(200,60,0.9,3,300000,200000,nothing,(0.5,2.5),[f1,f2])
end
calculate(s1)

=#
#=
vecsolve(s1)
v1,v2 = append!([f1,f2],vecsolve(s1)) |> shear 
m1 = moment(v1,3)
m2 = moment(v2,3)
x = collect(0:0.01:3)

plot(v1,(0,3))

yv1 = map(z->v1(z,3),x)
yv2 = map(z->v2(z,3),x)
ym1 = map(z->m1(z),x)
ym2 = map(z->m2(z),x)
#@show v1

#m1 = moment(v1,3)
#@show m1(2)
#plot(a,b) |> display
plot(x,ym1)
plot(x,ym2)
=#