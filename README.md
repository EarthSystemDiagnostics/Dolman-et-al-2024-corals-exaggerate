
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Dolman-et-al-2025-corals-exaggerate

This repository contains the data and code required to reproduce the
analyses in Dolman et al (2025) in three files:

- `01-processes-CH2K-gSST-data.qmd`
- `02-find-coral-pairs-attach-gSST.qmd`
- `03-Dolman-et-al-2024-corals-exaggerate-pairs.qmd`

The file `03-Dolman-et-al-2024-corals-exaggerate-pairs.qmd` reproduces
the main analysis and can be run using the pre-processed data in `/data`

To validate the full data processing steps, code can be run from the
following the two files:

`02-find-coral-pairs-attach-gSST.qmd` identifies the pairs of coral
records we used, does some interpolation across short gaps and attaches
pre-processed gridded SST data. This code can also be run using
pre-processed data in `/data`

`01-processes-CH2K-gSST-data.qmd` requires additional large raw data
files to be downloaded from their sources

The R serialisation of CoralHydro2k

<https://doi.org/10.25921/yp94-v135>
<https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/noaa-coral-35453/html>

and gridded SST data:

OISSTv2.1 at monthly resolution and 0.25’ available here:
<https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html>

File approximately 2GB, named “sst.mon.mean.nc”

ERSST at monthly resolution available here:
<https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/v6/one.big.file/ersst.v6.1850-2024.nc>

Dolman, A., McPartland, M., Felis, T. & Laepple, T. Coral records
exaggerate past decadal tropical climate variability. Preprint at
<https://doi.org/10.21203/rs.3.rs-3924954/v1> (2025).
