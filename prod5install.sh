#!/bin/bash

mkdir tmp
cd tmp
wget http://cta.irap.omp.eu/ctools/_downloads/csadd2caldb.py
wget https://zenodo.org/record/5499840/files/cta-prod5-zenodo-fitsonly-v0.1.zip
python csadd2caldb.py debug=yes
cd ..
rm -rf tmp
