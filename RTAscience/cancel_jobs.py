import os

runs = 20
start = 171

for i in range(runs):
    os.system(f'scancel {start+i}')