#!/usr/bin/env python

import sys
import math
n = float(sys.argv[1])
result = math.log(n)
print "Parsed args",':'.join(sys.argv)
print "log(%s) = %s" % (n,result)
