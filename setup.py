#!/usr/bin/python
# -*- coding: latin-1 -*-
from setuptools import setup, find_packages
setup( name='RTAscience',
       version='0.0.0',
       author='Ambra Di Piano',
       author_email='ambra.dipiano@inaf.it',
       packages=find_packages(),
       package_dir={ 'RTAscience': 'RTAscience' },
       include_package_data=True,
       license='BSD-3-Clause'
     )