# simple test of animate tool
# requires files/bucky*png files

a = animate("files/bucky*gif")
a.play()
a.delay(0.1)
a.back()
a.last()
a.first()
a.next()
a.previous()
a.frame(1)

print "all done ... type CTRL-D to exit Pizza.py"
