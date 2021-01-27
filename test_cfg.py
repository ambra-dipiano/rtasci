import yaml
import sys
import os
from os.path import join

cfgfile = sys.argv[1]

pypath = str(os.path.dirname(os.path.abspath(__file__)))  
configuration = open(join(pypath, cfgfile) )
cfg = yaml.load(configuration, Loader=yaml.FullLoader)

import simGRBpreparation

print(cfg)
