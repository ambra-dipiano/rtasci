import os

runs = 20
start = 491

for i in range(runs):
    os.system(f'scancel {start+i}')