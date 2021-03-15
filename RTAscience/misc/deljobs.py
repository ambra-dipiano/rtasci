import os
import sys

start = int(sys.argv[1])

runs = 20
for i in range(runs):
    os.system(f'scancel {start+i}')