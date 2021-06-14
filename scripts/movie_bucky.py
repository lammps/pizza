# movie of bucky-ball data

d = dump("dump.bucky")
d.set("$center = ((int($id)-1)/61+1)*61")
d.owrap("center")

r = raster(d)
r.file = "bucky"
r.acol([1,2,3],["blue","red","green"])
r.arad([1,2,3],[0.5,0.5,2.5])
r.box(1)

# spin

d.tselect.test("$t == 0")
r.pan(60,150,1,60,30,1)
r.all(0,25,0)

# slice

r.pan()
r.rotate(60,30)
r.select = "$x < -17 + (1-%g) * 34 and $type > 1"
r.all(0,25,25)

# timestep animation

d.tselect.all()
d.aselect.test("$type > 1")
r.rotate(60,30)
r.select = ""
r.all(50)
