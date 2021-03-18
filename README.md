### Environment **

To create a virtual environment with all required dependencies:

```bash
conda env create --name <envname> --file=environment.yaml
```

Note that you should already have anaconda installed: https://www.anaconda.com/

### Calibration database

To complete the environment be sure to download and install the correct IRFs (only prod2 comes with ctools installation). Public ones can be found here: https://www.cta-observatory.org/science/cta-performance/


### Configuration file

Under cfg you can find an example of configuration file. Description of each parameter is commented within. This file will serve as input when running the code. 

*__Required config parameters__*
Those with * are required only if offset = gw!

1. prepareGRBcatalog.py
    - setup: runid, scalefluxfactor
    - simulation: /
    - analysis: /
    - options: extract_data, set_ebl
    - path: data, catalog, ebl, model
2. simGRBcatalog.py
    - setup: simtype, runid, trials, start_count
    - simulation: caldb, irf, tobs, onset, delay, emin, emax, roi, offset
    - analysis: /
    - options: /
    - path: data, catalog, model, merger*
3. simWobble.py
    - setup: simtype, runid, trials, start_count
    - simulation: caldb, irf, tobs, emin, emax, roi, offset, nruns
    - analysis: /
    - options: /
    - path: data, catalog, model
4. simBkg.py
    - setup: simtype, trials, start_count
    - simulation: caldb, irf, tobs, emin, emax, roi, offset
    - analysis: /
    - options: /
    - path: data, catalog, model, bkg
5. pipelines/ctool3d_blindfit.py
    - setup: simtype, runid, trials, start_count
    - simulation: caldb, irf, tobs, onset, delay, emin, emax, roi, offset
    - analysis: tool, type, blind, binned, exposure, skypix, skyroifrac, smooth, maxsrc, sgmthresh, usepnt
    - options: plotsky
    - path: data, catalog, model, merger*
6. rtapipe.py
    - It combines one or more of the previous scripts, hence the required configuration will depend on their specifics.


    

### CALDB degradation
Be sure to have your calibration database installed under $CTOOLS/share. You can pass a single CALDB or a list, the degraded version will placed along side the nominal one. It will have the same suffix, replacing "prod" with "degr" (i.e., prod2 --> degr2).

```bash
python degradation_caldb.py --caldb prod3b pro3b-v2
```

Note: currently the code simply halves the affective area (and consequently renormalise the background rates).

### RTAscience
To extract spectra and lightcurves from the templates (one, a list or the entire sample):

```bash
python prepareGRBcatalog.py -f cfg/config.yaml
```

To run the simulation:

```bash
python simGRBcatalog.py -f cfg/config.yaml
```

To perform the analysis:

```bash
python pipeline/pipeline_name.py -f cfg/config.yaml
```
You are required to substitute pipeline_name.py with the chosen script (currently only one pipeline is available but there will be more in the future). 

All steps above may be run together with the following:

```bash
python rtapipe.py -f cfg/config.yaml
```
This will first extract all data specified in the configuration file, then simulate the entire sample, finally it will analyse each simulation. Files will be removed (by default) after the simulation so be sure to have enough space to store them. It will avoid running multiple simulation of the same sample if, i.e., you want to perform different types of analysis on it.

For more options on the above you can type:

```bash
python scriptName.py -h
```


<HR>
[**] subsceptible to changes 