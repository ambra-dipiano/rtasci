from setuptools import setup, find_packages

entry_points = {
	'console_scripts': [
		'simGRBcatalogWithRandomization = rtasci.simGRBcatalogWithRandomization:main',
		'simGRBcatalog = rtasci.simGRBcatalog:main'
     ]
}

setup(name='rtasci',
     author='Ambra Di Piano',
     author_email='ambra.dipiano@inaf.it',
     package_dir={'rtasci': 'rtasci'},
    packages=find_packages(),
    include_package_data=True,
    license='BSD-3-Clause',
)
