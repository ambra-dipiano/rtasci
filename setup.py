from setuptools import setup, find_packages

setup( 
     name='RTAscience',
     version='0.1.0',
     author='Ambra Di Piano',
     author_email='ambra.dipiano@inaf.it',
     packages=find_packages(),
     package_dir={ 'RTAscience': 'RTAscience' },
     include_package_data=True,
     license='BSD-3-Clause',
     python_requires=">=3.7"
)