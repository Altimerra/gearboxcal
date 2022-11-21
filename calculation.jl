begin
import Pkg
Pkg.activate(".")
include("calculator.jl")
import .Rolling, .Gears, .Shafts, .Bearings, .Keys
import PrettyTables
import Plots
end

begin
baseRollerRPM = 1.8
sheet = Rolling.Roller(30e-3, 6e-3, 1e-3, 40e-3,baseRollerRPM)
wire = Rolling.Roller(30e-3, 8e-3, 1e-3, 8e-3,baseRollerRPM*5.8)
brass = Rolling.Material(895e6, 0.49)

torque_sheet = Rolling.torque(sheet,brass)
power_sheet = Rolling.power(sheet,brass)

torque_wire = Rolling.torque(wire,brass)
power_wire = Rolling.power(wire,brass)
torque_sheet/torque_wire 
end

begin
power_motor = 60
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
end

begin
# Gear set 1
ga2 = gear(30e-3,12,2.413,true)
Gears.calculate(ga2)

gb2 = gear(50e-3,12,2.413,false)
Gears.calculate(ga2,gb2)

sd1 = Gears.centre_distance(ga2)

dia3 = 2*sd1 /(1+2)
gb3 = gear(dia3,12*2.075,2,true)
gb3.mod = 2.4
Gears.calculate2(gb3)
ga3 = gear(50e-3,12,2.413,false)
Gears.calculate(gb3,ga3)

sd11 = Gears.centre_distance(ga3)

end
# Gear set 2

gc1 = gear(30e-3,12,2.42,true)
Gears.calculate(gc1)

gc5 = gear(42.5e-3,12,2.5,true)
gc5.mod = 2.4
Gears.calculate2(gc5)

gb4 = gear(45e-3,4.826,2.413,true)
Gears.calculate(gb4)
gc4 = gear(81e-3,gb4.N,1,false)
Gears.calculate(gb4,gc4)

sd2 = Gears.centre_distance(gc4)

dia4 = 2*sd2 /(1+2.4)
sb22 = Gears.centre_distance(gc5)

gc5 = gear(58e-3,12,2.4,true)
gc5.mod = 2.9
Gears.calculate2(gc5)

gb5 = gear(60e-3,4.8,2.4,false)
Gears.calculate(gc5,gb5)

sd22 = Gears.centre_distance(gc5)

gc1 = gear(53e-3,4.8,2.45,true)
gc1.mod=2.9
Gears.calculate2(gc1)


gr1 = gear(40e-3,12,2.42,true)
gr1.mod = 2.9
Gears.calculate2(gr1)


macro saveplots(plots)
    return quote
        dir = "figures/" * $(string(plots))
        mkpath(dir)
        path1 = dir * "/shearx.png"
        path2 = dir * "/sheary.png"
        path3 = dir * "/momentx.png"
        path4 = dir * "/momenty.png"
        path5 = dir * "/all.png"
        Plots.savefig($plots[1], path1)        
        Plots.savefig($plots[2], path2)        
        Plots.savefig($plots[3], path3)        
        Plots.savefig($plots[4], path4)        
        Plots.savefig($plots[5], path5)        
    end
end


begin
shaft_len = 110e-3
sig = 700e6*0.36
tau = 0.5*sig
ka = 0.0
kb = 0.0

kc = 0.0
# Define the forces acting on shaft
fa2 = Shafts.Force(Gears.wt(ga2),Gears.wt(ga2)*tan(pres_angle),shaft_len/3)

sa = Shafts.Shaft(power_motor,ga2.N,ka,shaft_len/2,tau,sig,nothing,(0,shaft_len/2),[fa2])
sa_plots = Shafts.calculate(sa)
@saveplots sa_plots

fb4 = Shafts.Force(Gears.wt(gb4),Gears.wt(gb4)*tan(pres_angle),2*shaft_len/3)

sb = Shafts.Shaft(power_motor,gb5.N,kb,shaft_len,tau,sig,nothing,(0,shaft_len),[fa2,fb4])
sb_plots = Shafts.calculate(sb)
@saveplots sb_plots

sc = Shafts.Shaft(power_motor,gc5.N,kc,shaft_len,tau,sig,nothing,(0,shaft_len),[fb4])
sc_plots = Shafts.calculate(sc)
@saveplots sc_plots

rey=5526*( 1*sin(1.198) - tan(deg2rad(20))*cos(1.198) ) +  1989*( 1*cos(0.501) - tan(deg2rad(20))*sin(0.501) )
rex=1989*( 1*sin(1.198) - tan(deg2rad(20))*cos(1.198) ) +  5526*( 1*cos(0.501) - tan(deg2rad(20))*sin(0.501) )
fr1 = Shafts.Force(rex,rey,shaft_len/10)

sr = Shafts.Shaft(power_motor,gc5.N*54/42,kc,shaft_len*2/10,tau,sig,nothing,(0,shaft_len*2/10),[fr1])
sr_plots = Shafts.calculate(sr)
@saveplots sr_plots
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

sr.dia = ceil(findmax(sr.dia)[1];digits=3)
fr1 = sr.forces[lastindex(sc.forces)-1]
fr2 = sr.forces[lastindex(sc.forces)]


# Using same stress values from shaft
#convert diameter to mm
key(s) = Keys.Key(0,s.dia,0,0,sig,tau,Shafts.T(s),0)

ka = key(sa)
Keys.calculate(ka)
kb = key(sb)
Keys.calculate(kb)
kc = key(sc)
Keys.calculate(kc)
kr = key(sr)
Keys.calculate(kr)

# Weakening of shaft has to be considered

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
bearingCalc("r","1")
bearingCalc("r","2")

function tabler(type,datafunc,things,fields,headers,label)
    "$(type)_data = [ $datafunc(x) for x in $things ]" |> Meta.parse |> eval
    "$(type)_data = convert(Vector{Any}, $(type)_data)" |> Meta.parse |> eval  
    "pushfirst!($(type)_data,$fields)" |> Meta.parse |> eval
    "$(type)_data = reduce(hcat, $(type)_data)" |> Meta.parse |> eval
    "table_$(type) = PrettyTables.pretty_table(String, $(type)_data;
        header=$headers,backend = Val(:latex),hlines=:none,label=$label)" |> Meta.parse |> eval
    "open(\"tables/$type.tex\", \"w\") do file; write(file, table_$(type) );
        end" |> Meta.parse |> eval
end

function tabler2(type,datafunc,things,fields,headers)
    "$(type)_data = [ $datafunc(x) for x in $things ]" |> Meta.parse |> eval
    "$(type)_data = convert(Vector{Any}, $(type)_data)" |> Meta.parse |> eval  
    "pushfirst!($(type)_data,$fields)" |> Meta.parse |> eval
    "$(type)_data = reduce(hcat, $(type)_data)'" |> Meta.parse |> eval
    "table_$(type) = PrettyTables.pretty_table(String, $(type)_data;
        header=$headers,backend = Val(:latex),hlines=:none)" |> Meta.parse |> eval
    "open(\"tables/$type.tex\", \"w\") do file; write(file, table_$(type) );
        end" |> Meta.parse |> eval
end

data_gears(g) = [g.T ,g.mod, g.dia, g.sig, g.sige, g.siges, g.b[1]]
gears = [ga2,gb2,ga3,gb3,gb4,gc4,gb5,gc5,gc1,gr1]
label_gear = "Parameters of Gears"
fields_gears = ["Teeth","Module(mm)","Diameter(m)","Allowable stress(MPa)","Flextural Endurance(MPa)","Surface Endurance(MPa)","Face width(mm)"]
headers_gears=["Gear No.","Ga2","Gb2","Ga3","Gb3","Gb4","Gc4","Gb5","Gc5","Gc1","Gr1"]
tabler("gears","data_gears","gears","fields_gears","headers_gears", "label_gear")

data_shafts(s) = [s.dia, s.L, s.k]
shafts = [sa, sb, sc, sr]
label_shafts = "Parameters of Shafts"
fields_shafts = ["Diameter(m)", "Length(m)", "Diameter ratio"]
headers_shafts = ["","Sa", "Sb", "Sc", "Sr"]
tabler("shafts","data_shafts","shafts","fields_shafts","headers_shafts","label_shafts")

data_bearings(b) = [b.d, b.c, b.wr]
bearings = [ba1,ba2,bb1,bb2,bc1,bc2,br1,br2]
label_bear = "Parameters of Bearings"
fields_bearings = ["Min Diameter(m)", "Dynamic load rating", "Radial load(N)"]
headers_bearings = ["","Ba1","Ba2", "Bb1","Bb2", "Bc1","Bc2","Br1","Br2"]
tabler("bearings","data_bearings","bearings","fields_bearings","headers_bearings","label_bear")

data_keys(k) = [k.l, k.w, k.t, k.e]
keys = [ka,kb,kc,kr]
label_key = "Parameters of Keys"
fields_key = ["Length(m)", "Width(m)", "Thickness(m)","Shaft weakening factor"]
headers_key= ["","Ka","Kb","Kc","Kr"]
tabler("keys","data_keys","keys","fields_key","headers_key","label_key")


# Splines
splinestress(T,D) = 16*T/(pi*D^3)
splinestress(115,0.026)*1e-6
tau

splinedim(D) = (0.250*D,0.075*D,0.85*D)
splinedim(0.025+0.010)

(Gears.wt(gc1),Gears.wt(gb5))


rey=5526*( 1*sin(1.198) - tan(deg2rad(20))*cos(1.198) ) +  1989*( 1*cos(0.501) - tan(deg2rad(20))*sin(0.501) )
rex=1989*( 1*sin(1.198) - tan(deg2rad(20))*cos(1.198) ) +  5526*( 1*cos(0.501) - tan(deg2rad(20))*sin(0.501) )
res = rey^2+rex^2 |> sqrt
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
