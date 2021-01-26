import os
from lib.RTAVisualise import *
from os.path import join, isdir, isfile
from os import listdir

pypath = str(os.path.dirname(os.path.abspath(__file__)))  
datapath = pypath.replace('cta-sag-sci', 'DATA/outputs/crab')
pngpath = join(datapath, 'png')
if not isdir(pngpath):
    print('Creating png folder...')
    os.mkdir(pngpath)

tables = [f for f in listdir(datapath) if '.csv' in f]
print(len(tables))

for table in tables:
    data = pd.read_csv(join(datapath, table), sep=' ', )