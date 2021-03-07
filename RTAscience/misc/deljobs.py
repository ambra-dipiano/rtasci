import os
import sys

start = float(sys.argv[1])

runs = 20
for i in range(runs):
    os.system(f'scancel {start+i}')