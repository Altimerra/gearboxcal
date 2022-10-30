begin
import Pkg
Pkg.activate(".")
include("calculator.jl")
import .Rolling, .Gears, .Shafts, .Bearings
import PrettyTables
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
# TODO: Find real values
    (50,470), # static strength
    (80,700), # Flexural endurance limit
    (600,1050), # Surface endurace
    200000, # Elastic modulus 
    "mach", # Manufacturing
    nothing, # face width
    0 # Min face width factor
)


ga2 = gear(24e-3,60,5,true)
Gears.calculate(ga2)

gb2 = gear(50e-3,60,5,false)
Gears.calculate(ga2,gb2)


gb5 = gear(34e-3,60/5,6.67,true)
Gears.calculate(gb5)

gc5 =  gear(20e-3,60/5,6.67,false)
Gears.calculate(gb5,gc5)

sd1 = Gears.centre_distance(ga2)

Gears.centre_distance(gc5)

ga3 = gear(sd1,ga2.N,1,true)
Gears.calculate2(ga3)
gb3 = gear(sd1,ga3.N,1,false)
Gears.calculate(ga3,gb3)


sd2 = Gears.centre_distance(gb5)

gb4 = gear(sd2,gb5.N,1,true)
Gears.calculate2(gb4)

gc4 = gear(sd2,gb4.N,1,false)
Gears.calculate(gb4,gc4)
end

begin
shaft_len = 50e-3
sig = 320e6*0.6
tau = 0.3*sig
k = 0.7
# Define the forces acting on shaft
fa2 = Shafts.Force(Gears.wt(ga2),Gears.wt(ga2)*tan(pres_angle),shaft_len/3)

sa = Shafts.Shaft(power_motor,ga2.N,k,shaft_len,tau,sig,nothing,(0,shaft_len),[fa2])
Shafts.calculate(sa)

fb5 = Shafts.Force(Gears.wt(gb5),Gears.wt(gb5)*tan(pres_angle),2*shaft_len/3)

sb = Shafts.Shaft(power_motor,gb5.N,k,shaft_len,tau,sig,nothing,(0,shaft_len),[fa2,fb5])
Shafts.calculate(sb)

sc = Shafts.Shaft(power_motor,gc5.N,k,shaft_len,tau,sig,nothing,(0,shaft_len),[fb5])
Shafts.calculate(sc)
end

sa.dia = ceil(findmax(sa.dia)[1];digits=3)
fa1 = sa.forces[lastindex(sa.forces)-1]
fa2 = sa.forces[lastindex(sa.forces)]

sb.dia = ceil(findmax(sb.dia)[1];digits=3)
fb1 = sb.forces[lastindex(sb.forces)-1]
fb2 = sb.forces[lastindex(sb.forces)]

sc.dia = ceil(findmax(sc.dia)[1];digits=3)
fc1 = sc.forces[lastindex(sc.forces)-1]
fc2 = sc.forces[lastindex(sc.forces)]

bearing_life = 20*300*2 # hours

function bearingCalc(d,n)
Meta.parse("b$d$n = Bearings.Bearing(Shafts.resultant(f$d$n),0,s$d.N,bearing_life,s$d.dia,nothing,nothing)") |> eval
"Bearings.calculate(b$d$n)" |> Meta.parse |> eval
end

bearingCalc("a","1")
bearingCalc("a","2")
bearingCalc("b","1")
bearingCalc("b","2")
bearingCalc("c","1")
bearingCalc("c","2")

function tabler(type,datafunc,things,fields,headers)
    "$(type)_data = [ $datafunc(x) for x in $things ]" |> Meta.parse |> eval
    "$(type)_data = convert(Vector{Any}, $(type)_data)" |> Meta.parse |> eval  
    "pushfirst!($(type)_data,$fields)" |> Meta.parse |> eval
    "$(type)_data = reduce(hcat, $(type)_data)" |> Meta.parse |> eval
    "table_$(type) = PrettyTables.pretty_table(String, $(type)_data; header=$headers,backend = Val(:latex),hlines=:none)" |> Meta.parse |> eval
    "open(\"tables/$type.tex\", \"w\") do file; write(file, table_$(type) ); end" |> Meta.parse |> eval
end

function tabler2(type,datafunc,things,fields,headers)
    "$(type)_data = [ $datafunc(x) for x in $things ]" |> Meta.parse |> eval
    "$(type)_data = convert(Vector{Any}, $(type)_data)" |> Meta.parse |> eval  
    "pushfirst!($(type)_data,$fields)" |> Meta.parse |> eval
    "$(type)_data = reduce(hcat, $(type)_data)'" |> Meta.parse |> eval
    "table_$(type) = PrettyTables.pretty_table(String, $(type)_data; header=$headers,backend = Val(:latex),hlines=:none)" |> Meta.parse |> eval
    "open(\"tables/$type.tex\", \"w\") do file; write(file, table_$(type) ); end" |> Meta.parse |> eval
end

data_gears(g) = [g.T ,g.mod, g.dia, g.sig, g.sige, g.siges, g.b[1]]
gears = [ga2,gb2,ga3,gb3,gb4,gc4,gb4,gc4,gb5,gc5]
fields_gears = ["Teeth","Module(mm)","Diameter(m)","Allowable stress(MPa)","Flextural Endurance(MPa)","Surface Endurance(MPa)","Face width(mm)"]
headers_gears=["Gear No.","Ga2","Gb2","Ga3","Gb3","Gb4","Gc4","Gb4","Gc4","Gb5","Gc5"]
tabler("gears","data_gears","gears","fields_gears","headers_gears")

data_shafts(s) = [s.dia, s.L, s.k]
shafts = [sa, sb, sc]
fields_shafts = ["Diameter(m)", "Length(m)", "Diameter ratio"]
headers_shafts = ["","Sa", "Sb", "Sc"]
tabler("shafts","data_shafts","shafts","fields_shafts","headers_shafts")

data_bearings(b) = [b.d, b.c, b.wr]
bearings = [ba1,ba2,bb1,bb2,bc1,bc2]
fields_bearings = ["Diameter(m)", "Dynamic load rating", "Radial load(N)"]
headers_bearings = ["","Ba1","Ba2", "Bb1","Bb2", "Bc1","Bc2"]
tabler("bearings","data_bearings","bearings","fields_bearings","headers_bearings")



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