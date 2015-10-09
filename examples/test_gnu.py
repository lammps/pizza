# simple test of gnu tool
# creates tmp.eps

g = gnu()

g("plot sin(x) with lines")

a = range(10)
b = [3,6,2,5,7,3,6,5,3,1]

g.plot(a)
g.plot(a,b)
g.plot(a,b,b,a)
g.mplot(0,10,2,"tmp",a,b)

g.export("tmp.gnu",a,b)

g.select(2)
g.plot(a,b,b,a)
g.hide(1)

g.aspect(1.0)
g.title("My title","x-axis","y-axis")
g.xrange(1,8)
g.yrange(2,7)
g.label(5,4,"this is a test label")
g.curve(1,'g')
g.curve(2,'m')
g.ylog()

g.save("tmp")

print "all done ... type CTRL-D to exit Pizza.py"
