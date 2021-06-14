# movie of Potts model grain growth

m = mdump("tmp.grain")
m.map(2,"spin")
m.etype = "spin"

r = raster(m)
r.rotate(0,0)
r.zoom(1.5)
r.all()
