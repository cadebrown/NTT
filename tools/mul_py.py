#!/usr/bin/env python3

import sys
import time

A = int(open(sys.argv[1]).read() if "." in sys.argv[1] else sys.argv[1], 16)
B = int(open(sys.argv[2]).read() if "." in sys.argv[2] else sys.argv[2], 16)

st = time.time()
C = A * B
st = time.time() - st

print ("time: %.3f" % (st, ), file=sys.stderr)
print (hex(C))


