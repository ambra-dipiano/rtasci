# Changelog

## **v.0.1.0**
- script to degrade caldb
- script to extract data GRB afterglow template, normalise the template, apply EBL when required
- scripts to simulate: grb afterglow, crab wobble, empty fields
- ctools pipeline: blind-search + full-FOV unbinned analysis
- simulation scripts now sort events in TIME and reindex them

#### *Known issues*
- the grb afterglow simulation (simGRBcatalog.py) for serendipitous discovery (onset != 0) is bugged. 