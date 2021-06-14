# create series of pressure/volume plots for flow around obstacle
# montage them with flow images to make flow/plot movie

lg = log("log.flow")
s,p,v = lg.get("Step","Press","Vol")

g = gnu()
g.aspect(1.0)
g.xrange(0,25000)
g.yrange(0,3)
g.title("Pressure","TimeSteps","Pressure")
g.mplot(0,251,1,"plotp",s,p)
g.yrange(700,950)
g.curve(1,'g')
g.title("Volume","TimeSteps","Volume")
g.mplot(0,251,1,"plotv",s,v)

i = image()
i.convert("plotp*eps","plotp*png","-crop 512x512+0+225")
i.convert("plotv*eps","plotv*png","-crop 512x512+0+225")
i.montage("-geometry 256x256","flow*png","plotp*png","plotv*png",
          "flow_plot*png")
i.convert("flow_plot*png","flow_plot.mpg")
