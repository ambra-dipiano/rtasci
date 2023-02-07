# *****************************************************************************
# Copyright (C) 2020 INAF
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *****************************************************************************

from setuptools import setup, find_packages

setup(name='RTAscience',
     author='Ambra Di Piano',
     author_email='ambra.dipiano@inaf.it',
     package_dir={'RTAscience': 'RTAscience'},
    packages=find_packages(),
    include_package_data=True,
    license='BSD-3-Clause',
)
