#!/usr/bin/env python

def f(x):
    import math
    return 10*math.e**(math.log(0.5)/5.27 * x)
def radiationExposure (start,stop,step):
    height=0
    x=start
    while (x<stop):
	height=f(x)+height
	x=x+step
    return height*step
print radiationExposure(40,100,1.5)
