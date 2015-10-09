# simple test of matlab tool
# creates tmp.eps

m = matlab()

a = range(10)
b = [3,6,2,5,7,3,6,5,3,1]

m.plot(a)
m.plot(a,b)
m.plot(a,b,b,a)
m.mplot(0,10,2,"tmp",a,b)

m("3*400 - 50")

m.export("tmp.gnu",a,b)

m.select(2)
m.plot(a,b,b,a)
m.hide(1)

m.aspect(1.0)
m.title("My title","x-axis","y-axis")
m.xrange(1,8)
m.yrange(2,7)
m.label(5,4,"this is a test label")
m.curve(1,'g')
m.curve(2,'m')
#m.ylog()

m.save("tmp")

print "all done ... type CTRL-D to exit Pizza.py"
