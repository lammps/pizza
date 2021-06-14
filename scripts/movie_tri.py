# movie of triangular particle data

d = dump("dump.tri")
d.set("$center = ((int($id)-1)/36+1)*36")
d.owrap("center")

r = raster(d)
r.bg("white")
r.file = "tri"
r.acol([1,2],["blue","red"])
r.arad([1,2],0.5)
r.box(1)

# spin

d.tselect.test("$t == 0")
r.pan(60,150,1,60,30,1)
r.all(0,25,0)

# slice

r.pan()
r.rotate(60,30)
r.select = "$x < -12 + (1-%g) * 24"
r.all(0,25,25)

# timestep animation

d.tselect.all()
r.rotate(60,30)
r.select = ""
r.all(50)
