# movie of melting LJ solid

a = dump("dump.melt")
a.tselect.test("$t == 0")
a.scale()
a.set("$ix = int($x * 4)")
a.set("$iy = int($y * 4)")
a.set("$iz = int($z * 4)")
a.set("$type = ($ix + $iy + $iz) % 2 + 1")
a.unscale()
a.tselect.all()
a.clone(0,"type")

r = raster(a)
r.acol([1,2],["red","green"])
r.arad([1,2],0.5)
r.file = "melt"
r.pan(130,25,1,60,135,0.6)

r.all()
