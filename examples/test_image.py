# simple test of image tool
# requires files/bucky*.png

i = image("files/bucky*.gif")
i.convert("files/bucky*.gif","tmp*.png")
i.montage("","files/bucky*.gif","tmp*.png","tmpnew*.gif")
i.view("*.gif")

print "all done ... type CTRL-D to exit Pizza.py"
