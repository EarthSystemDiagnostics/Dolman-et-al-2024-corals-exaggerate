
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Dolman-et-al-2026-corals-exaggerate

This repository contains the data and code required to reproduce the
analyses in Dolman et al (2026) in three files:

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

## Data summary

A summary of the records from CoralHydro2k that were used in the main
analysis of Dolman et al. (2026)

``` r
table_summary_of_paired_records <- readRDS("data/table_summary_of_paired_records.RDS")
knitr::kable(table_summary_of_paired_records,
             caption = "")
```

| paleoData_variableName | paleoData_TSid | pair | dist_km | paleoData_archiveSpecies | duration | first_joint_year | last_joint_year | geo_ocean | geo_siteName | pub1_firstauthor | pub1_year |
|:---|:---|:---|---:|:---|---:|---:|---:|:---|:---|:---|---:|
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-FL17DTO02_SrCa | 0.16 | Siderastrea siderea | 131 | 1837 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO02_SrCa | DE14DTO01_SrCa-FL17DTO02_SrCa | 0.16 | Siderastrea siderea | 131 | 1837 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-DE14DTO02_SrCa | 11.48 | Siderastrea siderea | 125 | 1779 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO02_SrCa | DE14DTO01_SrCa-DE14DTO02_SrCa | 11.48 | Siderastrea siderea | 125 | 1779 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-FL17DTO02_SrCa | 0.06 | Orbicella faveolata | 116 | 1893 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL17DTO02_SrCa | FL17DTO01_SrCa-FL17DTO02_SrCa | 0.06 | Siderastrea siderea | 116 | 1893 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL17DTO01_SrCa | 0.57 | Siderastrea siderea | 112 | 1896 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO01_SrCa | DE14DTO04_SrCa-FL17DTO01_SrCa | 0.57 | Orbicella faveolata | 112 | 1896 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL17DTO02_SrCa | 0.63 | Siderastrea siderea | 112 | 1896 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO02_SrCa | DE14DTO04_SrCa-FL17DTO02_SrCa | 0.63 | Siderastrea siderea | 112 | 1896 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | MU18GSI01_SrCa | MU18GSI01_SrCa-MU18NPI01_SrCa | 38.90 | Porites sp. | 108 | 1908 | 2016 | Pacific Ocean | Gili Selang, Bali, Indonesia | Murty, S. A. | 2018 |
| SrCa | MU18NPI01_SrCa | MU18GSI01_SrCa-MU18NPI01_SrCa | 38.90 | Porites sp. | 108 | 1908 | 2016 | Indian Ocean | Nusa Penida, Indonesia | Murty, S. A. | 2018 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-MA08DTO01_SrCa | 11.32 | Siderastrea siderea | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | MA08DTO01_SrCa | DE14DTO04_SrCa-MA08DTO01_SrCa | 11.32 | Siderastrea siderea | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-MA08DTO01_SrCa | 11.39 | Orbicella faveolata | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | MA08DTO01_SrCa | FL17DTO01_SrCa-MA08DTO01_SrCa | 11.39 | Siderastrea siderea | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | FL17DTO02_SrCa | FL17DTO02_SrCa-MA08DTO01_SrCa | 11.40 | Siderastrea siderea | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | MA08DTO01_SrCa | FL17DTO02_SrCa-MA08DTO01_SrCa | 11.40 | Siderastrea siderea | 93 | 1900 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-FL17DTO01_SrCa | 0.12 | Siderastrea siderea | 75 | 1893 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO01_SrCa | DE14DTO01_SrCa-FL17DTO01_SrCa | 0.12 | Orbicella faveolata | 75 | 1893 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-DE14DTO04_SrCa | 0.51 | Siderastrea siderea | 72 | 1896 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO04_SrCa | DE14DTO01_SrCa-DE14DTO04_SrCa | 0.51 | Siderastrea siderea | 72 | 1896 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-MA08DTO01_SrCa | 11.49 | Siderastrea siderea | 68 | 1900 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | MA08DTO01_SrCa | DE14DTO01_SrCa-MA08DTO01_SrCa | 11.49 | Siderastrea siderea | 68 | 1900 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | DE14DTO02_SrCa | DE14DTO02_SrCa-FL17DTO02_SrCa | 11.40 | Siderastrea siderea | 67 | 1837 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO02_SrCa | DE14DTO02_SrCa-FL17DTO02_SrCa | 11.40 | Siderastrea siderea | 67 | 1837 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-DE14DTO04_SrCa | 0.51 | Siderastrea siderea | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO04_SrCa | DE14DTO03_SrCa-DE14DTO04_SrCa | 0.51 | Siderastrea siderea | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL17DTO01_SrCa | 0.12 | Siderastrea siderea | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO01_SrCa | DE14DTO03_SrCa-FL17DTO01_SrCa | 0.12 | Orbicella faveolata | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL17DTO02_SrCa | 0.16 | Siderastrea siderea | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO02_SrCa | DE14DTO03_SrCa-FL17DTO02_SrCa | 0.16 | Siderastrea siderea | 64 | 1945 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | ZI16ROD01_SrCa | ZI16ROD01_SrCa-ZI16ROD02_SrCa | 0.69 | Porites sp. | 60 | 1946 | 2006 | Indian Ocean | Rodrigues, Republic of Mauritius | Zinke, J. | 2016 |
| SrCa | ZI16ROD02_SrCa | ZI16ROD01_SrCa-ZI16ROD02_SrCa | 0.69 | Porites lobata | 60 | 1946 | 2006 | Indian Ocean | Rodrigues, Republic of Mauritius | Zinke, J. | 2016 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-MA08DTO01_SrCa | 11.49 | Siderastrea siderea | 48 | 1945 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | MA08DTO01_SrCa | DE14DTO03_SrCa-MA08DTO01_SrCa | 11.49 | Siderastrea siderea | 48 | 1945 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | GR13MAD01_SrCa | GR13MAD01_SrCa-GR13MAD02_SrCa | 0.74 | Porites lutea | 43 | 1964 | 2007 | Indian Ocean | Nosy Boraha, Madagascar (formerly Ile Sainte-Marie) | Grove, C. A. | 2013 |
| SrCa | GR13MAD02_SrCa | GR13MAD01_SrCa-GR13MAD02_SrCa | 0.74 | Porites lutea | 43 | 1964 | 2007 | Indian Ocean | Nosy Boraha, Madagascar (formerly Ile Sainte-Marie) | Grove, C. A. | 2013 |
| SrCa | SM06LKF01_SrCa | SM06LKF01_SrCa-SM06LKF02_SrCa | 0.00 | Orbicella faveolata | 42 | 1961 | 2003 | Atlantic Ocean | Looe Key, Florida, USA | Smith, J. M. | 2006 |
| SrCa | SM06LKF02_SrCa | SM06LKF01_SrCa-SM06LKF02_SrCa | 0.00 | Orbicella faveolata | 42 | 1961 | 2003 | Atlantic Ocean | Looe Key, Florida, USA | Smith, J. M. | 2006 |
| SrCa | BA04FIJ02_SrCa | BA04FIJ02_SrCa-LI04FIJ01_SrCa | 0.00 | Diploastrea heliopora | 35 | 1962 | 1997 | Pacific Ocean | Savusavu Bay, Vanua Levu, Fiji | Bagnato, S. | 2004 |
| SrCa | LI04FIJ01_SrCa | BA04FIJ02_SrCa-LI04FIJ01_SrCa | 0.00 | Porites lutea | 35 | 1962 | 1997 | Pacific Ocean | Vanua Levu, Fiji | Linsley, B. K. | 2004 |
| SrCa | JI18GAL01_SrCa | JI18GAL01_SrCa-JI18GAL02_SrCa | 0.00 | Porites lobata | 34 | 1976 | 2010 | Pacific Ocean | Shark Bay, Wolf Island, Ecuador | Jimenez, G. | 2018 |
| SrCa | JI18GAL02_SrCa | JI18GAL01_SrCa-JI18GAL02_SrCa | 0.00 | Porites lobata | 34 | 1976 | 2010 | Pacific Ocean | Shark Bay, Wolf Island, Ecuador | Jimenez, G. | 2018 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL18DTO02_SrCa | 0.53 | Siderastrea siderea | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO02_SrCa | DE14DTO03_SrCa-FL18DTO02_SrCa | 0.53 | Orbicella faveolata | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL18DTO02_SrCa | 0.04 | Siderastrea siderea | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO02_SrCa | DE14DTO04_SrCa-FL18DTO02_SrCa | 0.04 | Orbicella faveolata | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-FL18DTO02_SrCa | 0.58 | Orbicella faveolata | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO02_SrCa | FL17DTO01_SrCa-FL18DTO02_SrCa | 0.58 | Orbicella faveolata | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO02_SrCa | FL17DTO02_SrCa-FL18DTO02_SrCa | 0.64 | Siderastrea siderea | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO02_SrCa | FL17DTO02_SrCa-FL18DTO02_SrCa | 0.64 | Orbicella faveolata | 29 | 1980 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | BO14HTI01_SrCa | BO14HTI01_SrCa-BO14HTI02_SrCa | 0.00 | Porites lutea | 28 | 1977 | 2005 | Pacific Ocean | Hon Tre Island, Vietnam | Bolton, A. | 2014 |
| SrCa | BO14HTI02_SrCa | BO14HTI01_SrCa-BO14HTI02_SrCa | 0.00 | Porites lutea | 28 | 1977 | 2005 | Pacific Ocean | Hon Tre Island, Vietnam | Bolton, A. | 2014 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL18DTO01_SrCa | 0.44 | Siderastrea siderea | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO01_SrCa | DE14DTO03_SrCa-FL18DTO01_SrCa | 0.44 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL18DTO01_SrCa | 0.36 | Siderastrea siderea | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO01_SrCa | DE14DTO04_SrCa-FL18DTO01_SrCa | 0.36 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-FL18DTO01_SrCa | 0.41 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO01_SrCa | FL17DTO01_SrCa-FL18DTO01_SrCa | 0.41 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO02_SrCa | FL17DTO02_SrCa-FL18DTO01_SrCa | 0.46 | Siderastrea siderea | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO01_SrCa | FL17DTO02_SrCa-FL18DTO01_SrCa | 0.46 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO01_SrCa | FL18DTO01_SrCa-FL18DTO02_SrCa | 0.35 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO02_SrCa | FL18DTO01_SrCa-FL18DTO02_SrCa | 0.35 | Orbicella faveolata | 26 | 1983 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE14DTO01_SrCa | DE14DTO01_SrCa-DE14DTO03_SrCa | 0.00 | Siderastrea siderea | 23 | 1945 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | DE14DTO03_SrCa | DE14DTO01_SrCa-DE14DTO03_SrCa | 0.00 | Siderastrea siderea | 23 | 1945 | 1968 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO02_SrCa | FL17DTO02_SrCa-FL18DTO04_SrCa | 4.58 | Siderastrea siderea | 23 | 1989 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO04_SrCa | FL17DTO02_SrCa-FL18DTO04_SrCa | 4.58 | Orbicella faveolata | 23 | 1989 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | HE18COC01_SrCa | HE18COC01_SrCa-HE18COC02_SrCa | 1.01 | Porites sp. | 22 | 1988 | 2010 | Indian Ocean | Cocos (Keeling) Islands, Australia | Hennekam, R. | 2018 |
| SrCa | HE18COC02_SrCa | HE18COC01_SrCa-HE18COC02_SrCa | 1.01 | Porites sp. | 22 | 1988 | 2010 | Indian Ocean | Cocos (Keeling) Islands, Australia | Hennekam, R. | 2018 |
| SrCa | CH18YOA01_SrCa | CH18YOA01_SrCa-CH18YOA02_SrCa | 0.00 | Porites sp. | 19 | 1994 | 2013 | Pacific Ocean | Lingyang Reef, Yongle Atoll | Chen, T. | 2018 |
| SrCa | CH18YOA02_SrCa | CH18YOA01_SrCa-CH18YOA02_SrCa | 0.00 | Porites sp. | 19 | 1994 | 2013 | Pacific Ocean | Lingyang Reef, Yongle Atoll | Chen, T. | 2018 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL18DTO04_SrCa | 4.74 | Siderastrea siderea | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO04_SrCa | DE14DTO03_SrCa-FL18DTO04_SrCa | 4.74 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL18DTO04_SrCa | 5.07 | Siderastrea siderea | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO04_SrCa | DE14DTO04_SrCa-FL18DTO04_SrCa | 5.07 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-FL18DTO04_SrCa | 4.62 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO04_SrCa | FL17DTO01_SrCa-FL18DTO04_SrCa | 4.62 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO01_SrCa | FL18DTO01_SrCa-FL18DTO04_SrCa | 4.71 | Orbicella faveolata | 19 | 1989 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO04_SrCa | FL18DTO01_SrCa-FL18DTO04_SrCa | 4.71 | Orbicella faveolata | 19 | 1989 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO02_SrCa | FL18DTO02_SrCa-FL18DTO04_SrCa | 5.06 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO04_SrCa | FL18DTO02_SrCa-FL18DTO04_SrCa | 5.06 | Orbicella faveolata | 19 | 1989 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | NU11PAL01_SrCa | NU11PAL01_SrCa-SA19PAL01_SrCa | 1.57 | Porites lutea | 17 | 1981 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Nurhati, I. S. | 2011 |
| SrCa | SA19PAL01_SrCa | NU11PAL01_SrCa-SA19PAL01_SrCa | 1.57 | Porites sp. | 17 | 1981 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| SrCa | FL18DTO03_SrCa | FL18DTO03_SrCa-FL18DTO04_SrCa | 0.40 | Orbicella faveolata | 15 | 1998 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO04_SrCa | FL18DTO03_SrCa-FL18DTO04_SrCa | 0.40 | Orbicella faveolata | 15 | 1998 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL17DTO02_SrCa | FL17DTO02_SrCa-FL18DTO03_SrCa | 4.98 | Siderastrea siderea | 14 | 1998 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO03_SrCa | FL17DTO02_SrCa-FL18DTO03_SrCa | 4.98 | Orbicella faveolata | 14 | 1998 | 2012 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO02_SrCa | FL18DTO02_SrCa-MA08DTO01_SrCa | 11.28 | Orbicella faveolata | 13 | 1980 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | MA08DTO01_SrCa | FL18DTO02_SrCa-MA08DTO01_SrCa | 11.28 | Siderastrea siderea | 13 | 1980 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | NU11PAL01_SrCa | NU11PAL01_SrCa-SA19PAL02_SrCa | 1.57 | Porites lutea | 13 | 1985 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Nurhati, I. S. | 2011 |
| SrCa | SA19PAL02_SrCa | NU11PAL01_SrCa-SA19PAL02_SrCa | 1.57 | Porites sp. | 13 | 1985 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| SrCa | SA19PAL01_SrCa | SA19PAL01_SrCa-SA19PAL02_SrCa | 0.00 | Porites sp. | 13 | 1985 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| SrCa | SA19PAL02_SrCa | SA19PAL01_SrCa-SA19PAL02_SrCa | 0.00 | Porites sp. | 13 | 1985 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| SrCa | RE19GBR04_SrCa | RE19GBR04_SrCa-RE19GBR05_SrCa | 13.32 | Porites sp. | 12 | 1993 | 2005 | Pacific Ocean | Nomad Reef, Australia | Reed, E. V. | 2019 |
| SrCa | RE19GBR05_SrCa | RE19GBR04_SrCa-RE19GBR05_SrCa | 13.32 | Porites sp. | 12 | 1993 | 2005 | Pacific Ocean | Clerke Reef, Australia | Reed, E. V. | 2019 |
| SrCa | DE14DTO02_SrCa | DE14DTO02_SrCa-FL17DTO01_SrCa | 11.39 | Siderastrea siderea | 11 | 1893 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL17DTO01_SrCa | DE14DTO02_SrCa-FL17DTO01_SrCa | 11.39 | Orbicella faveolata | 11 | 1893 | 1904 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | DE14DTO03_SrCa | DE14DTO03_SrCa-FL18DTO03_SrCa | 5.14 | Siderastrea siderea | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO03_SrCa | DE14DTO03_SrCa-FL18DTO03_SrCa | 5.14 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE14DTO04_SrCa | DE14DTO04_SrCa-FL18DTO03_SrCa | 5.47 | Siderastrea siderea | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | DeLong, K. | 2014 |
| SrCa | FL18DTO03_SrCa | DE14DTO04_SrCa-FL18DTO03_SrCa | 5.47 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | DE16RED01_SrCa | DE16RED01_SrCa-MU18RED04_SrCa | 28.27 | Porites sp. | 11 | 1998 | 2009 | Red Sea | Red Sea | DeCarlo, T. M. | 2016 |
| SrCa | MU18RED04_SrCa | DE16RED01_SrCa-MU18RED04_SrCa | 28.27 | Porites lutea | 11 | 1998 | 2009 | Red Sea | Coral Gardens, Red Sea | Murty, S. A. | 2018 |
| SrCa | FL17DTO01_SrCa | FL17DTO01_SrCa-FL18DTO03_SrCa | 5.02 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2017 |
| SrCa | FL18DTO03_SrCa | FL17DTO01_SrCa-FL18DTO03_SrCa | 5.02 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO01_SrCa | FL18DTO01_SrCa-FL18DTO03_SrCa | 5.11 | Orbicella faveolata | 11 | 1998 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO03_SrCa | FL18DTO01_SrCa-FL18DTO03_SrCa | 5.11 | Orbicella faveolata | 11 | 1998 | 2008 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO01_SrCa | FL18DTO01_SrCa-MA08DTO01_SrCa | 11.09 | Orbicella faveolata | 11 | 1983 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | MA08DTO01_SrCa | FL18DTO01_SrCa-MA08DTO01_SrCa | 11.09 | Siderastrea siderea | 11 | 1983 | 1993 | Atlantic Ocean | Dry Tortugas, Florida, USA | Maupin, C | 2008 |
| SrCa | FL18DTO02_SrCa | FL18DTO02_SrCa-FL18DTO03_SrCa | 5.45 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | FL18DTO03_SrCa | FL18DTO02_SrCa-FL18DTO03_SrCa | 5.45 | Orbicella faveolata | 11 | 1998 | 2009 | Atlantic Ocean | Dry Tortugas, Florida, USA | Flannery, J. | 2018 |
| SrCa | RE19GBR01_SrCa | RE19GBR01_SrCa-RE19GBR02_SrCa | 26.34 | Porites sp. | 11 | 1979 | 1990 | Pacific Ocean | Eel Reef, Australia | Reed, E. V. | 2019 |
| SrCa | RE19GBR02_SrCa | RE19GBR01_SrCa-RE19GBR02_SrCa | 26.34 | Porites sp. | 11 | 1979 | 1990 | Pacific Ocean | Portland Roads, Australia | Reed, E. V. | 2019 |
| SrCa | KR20SAR01_SrCa | KR20SAR01_SrCa-KR20SAR02_SrCa | 9.82 | Porites sp. | 10 | 2006 | 2016 | Pacific Ocean | Sarawak, Malaysia | Krawczyk, H. | 2020 |
| SrCa | KR20SAR02_SrCa | KR20SAR01_SrCa-KR20SAR02_SrCa | 9.82 | Porites sp. | 10 | 2006 | 2016 | Pacific Ocean | Sarawak, Malaysia | Krawczyk, H. | 2020 |
| SrCa | XU15BVI01_SrCa | XU15BVI01_SrCa-XU15BVI03_SrCa | 0.00 | Pseudodiploria strigosa | 10 | 1998 | 2008 | Atlantic Ocean | Anegada, British Virgin Islands | Xu, Y. | 2015 |
| SrCa | XU15BVI03_SrCa | XU15BVI01_SrCa-XU15BVI03_SrCa | 0.00 | Pseudodiploria strigosa | 10 | 1998 | 2008 | Atlantic Ocean | Anegada, British Virgin Islands | Xu, Y. | 2015 |
| d18O | CH03LOM01_d18O | CH03LOM01_d18O-MU18NPI01_d18O | 46.46 | Porites sp. | 166 | 1824 | 1990 | Pacific Ocean | Padang Bai, Bali, Indonesia | Charles, C. D. | 2003 |
| d18O | MU18NPI01_d18O | CH03LOM01_d18O-MU18NPI01_d18O | 46.46 | Porites sp. | 166 | 1824 | 1990 | Indian Ocean | Nusa Penida, Indonesia | Murty, S. A. | 2018 |
| d18O | LI00RAR01_d18O | LI00RAR01_d18O-LI06RAR02_d18O | 0.33 | Porites lutea | 122 | 1875 | 1997 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2000 |
| d18O | LI06RAR02_d18O | LI00RAR01_d18O-LI06RAR02_d18O | 0.33 | Porites lutea | 122 | 1875 | 1997 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2004 |
| d18O | MU18GSI01_d18O | MU18GSI01_d18O-MU18NPI01_d18O | 38.90 | Porites sp. | 111 | 1905 | 2016 | Pacific Ocean | Gili Selang, Bali, Indonesia | Murty, S. A. | 2018 |
| d18O | MU18NPI01_d18O | MU18GSI01_d18O-MU18NPI01_d18O | 38.90 | Porites sp. | 111 | 1905 | 2016 | Indian Ocean | Nusa Penida, Indonesia | Murty, S. A. | 2018 |
| d18O | OS14RIP01_d18O | OS14RIP01_d18O-OS14UCP01_d18O | 14.83 | Porites lutea | 109 | 1899 | 2008 | Pacific Ocean | Rock Islands, Palau | Osborne, M. C. | 2014 |
| d18O | OS14UCP01_d18O | OS14RIP01_d18O-OS14UCP01_d18O | 14.83 | Porites lutea | 109 | 1899 | 2008 | Pacific Ocean | Ulong Channel, Palau | Osborne, M. C. | 2014 |
| d18O | CO93TAR01_d18O | CO93TAR01_d18O-UR00MAI01_d18O | 46.56 | Porites lobata | 96 | 1894 | 1990 | Pacific Ocean | Tarawa Atoll, Republic of Kiribati | Cole, J. E. | 1993 |
| d18O | UR00MAI01_d18O | CO93TAR01_d18O-UR00MAI01_d18O | 46.56 | Porites sp. | 96 | 1894 | 1990 | Pacific Ocean | Maiana, Republic of Kiribati | Urban, F. E. | 2000 |
| d18O | LI06RAR01_d18O | LI06RAR01_d18O-LI06RAR02_d18O | 0.00 | Porites lutea | 93 | 1907 | 2000 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2004 |
| d18O | LI06RAR02_d18O | LI06RAR01_d18O-LI06RAR02_d18O | 0.00 | Porites lutea | 93 | 1907 | 2000 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2004 |
| d18O | LI00RAR01_d18O | LI00RAR01_d18O-LI06RAR01_d18O | 0.33 | Porites lutea | 90 | 1907 | 1997 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2000 |
| d18O | LI06RAR01_d18O | LI00RAR01_d18O-LI06RAR01_d18O | 0.33 | Porites lutea | 90 | 1907 | 1997 | Pacific Ocean | Rarotonga, Cook Islands | Linsley, B. K. | 2004 |
| d18O | CH03LOM01_d18O | CH03LOM01_d18O-MU18GSI01_d18O | 27.24 | Porites sp. | 85 | 1905 | 1990 | Pacific Ocean | Padang Bai, Bali, Indonesia | Charles, C. D. | 2003 |
| d18O | MU18GSI01_d18O | CH03LOM01_d18O-MU18GSI01_d18O | 27.24 | Porites sp. | 85 | 1905 | 1990 | Pacific Ocean | Gili Selang, Bali, Indonesia | Murty, S. A. | 2018 |
| d18O | CA13DIA01_d18O | CA13DIA01_d18O-CA13PEL01_d18O | 50.77 | Orbicella faveolata | 66 | 1940 | 2007 | Atlantic Ocean | Diamond Caye, Utila, Honduras | Carilli, J. E. | 2013 |
| d18O | CA13PEL01_d18O | CA13DIA01_d18O-CA13PEL01_d18O | 50.77 | Orbicella faveolata | 66 | 1940 | 2007 | Atlantic Ocean | Cayos Cochinos, Honduras | Carilli, J. E. | 2013 |
| d18O | BA04FIJ02_d18O | BA04FIJ02_d18O-LI04FIJ01_d18O | 0.00 | Diploastrea heliopora | 57 | 1940 | 1997 | Pacific Ocean | Savusavu Bay, Vanua Levu, Fiji | Bagnato, S. | 2004 |
| d18O | LI04FIJ01_d18O | BA04FIJ02_d18O-LI04FIJ01_d18O | 0.00 | Porites lutea | 57 | 1940 | 1997 | Pacific Ocean | Vanua Levu, Fiji | Linsley, B. K. | 2004 |
| d18O | OS13NGP01_d18O | OS13NGP01_d18O-OS14RIP01_d18O | 16.04 | Porites lutea | 45 | 1963 | 2008 | Pacific Ocean | Ngaragabel, Palau | Osborne, M. C. | 2013 |
| d18O | OS14RIP01_d18O | OS13NGP01_d18O-OS14RIP01_d18O | 16.04 | Porites lutea | 45 | 1963 | 2008 | Pacific Ocean | Rock Islands, Palau | Osborne, M. C. | 2014 |
| d18O | OS13NGP01_d18O | OS13NGP01_d18O-OS14UCP01_d18O | 24.39 | Porites lutea | 45 | 1963 | 2008 | Pacific Ocean | Ngaragabel, Palau | Osborne, M. C. | 2013 |
| d18O | OS14UCP01_d18O | OS13NGP01_d18O-OS14UCP01_d18O | 24.39 | Porites lutea | 45 | 1963 | 2008 | Pacific Ocean | Ulong Channel, Palau | Osborne, M. C. | 2014 |
| d18O | CO03PAL04_d18O | CO03PAL04_d18O-CO03PAL05_d18O | 0.00 | Porites lutea | 43 | 1405 | 1448 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL05_d18O | CO03PAL04_d18O-CO03PAL05_d18O | 0.00 | Porites lutea | 43 | 1405 | 1448 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | SM06LKF01_d18O | SM06LKF01_d18O-SM06LKF02_d18O | 0.00 | Orbicella faveolata | 42 | 1961 | 2003 | Atlantic Ocean | Looe Key, Florida, USA | Smith, J. M. | 2006 |
| d18O | SM06LKF02_d18O | SM06LKF01_d18O-SM06LKF02_d18O | 0.00 | Orbicella faveolata | 42 | 1961 | 2003 | Atlantic Ocean | Looe Key, Florida, USA | Smith, J. M. | 2006 |
| d18O | CO03PAL04_d18O | CO03PAL04_d18O-CO03PAL06_d18O | 0.00 | Porites lutea | 32 | 1412 | 1443 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL06_d18O | CO03PAL04_d18O-CO03PAL06_d18O | 0.00 | Porites lutea | 32 | 1412 | 1443 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL05_d18O | CO03PAL05_d18O-CO03PAL06_d18O | 0.00 | Porites lutea | 32 | 1412 | 1443 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL06_d18O | CO03PAL05_d18O-CO03PAL06_d18O | 0.00 | Porites lutea | 32 | 1412 | 1443 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL03_d18O | CO03PAL03_d18O-CO03PAL10_d18O | 0.00 | Porites lutea | 31 | 1327 | 1357 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL10_d18O | CO03PAL03_d18O-CO03PAL10_d18O | 0.00 | Porites lutea | 31 | 1327 | 1357 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL08_d18O | CO03PAL08_d18O-CO03PAL09_d18O | 0.00 | Porites lutea | 31 | 1665 | 1696 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL09_d18O | CO03PAL08_d18O-CO03PAL09_d18O | 0.00 | Porites lutea | 31 | 1665 | 1696 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | BO14HTI01_d18O | BO14HTI01_d18O-BO14HTI02_d18O | 0.00 | Porites lutea | 28 | 1977 | 2005 | Pacific Ocean | Hon Tre Island, Vietnam | Bolton, A. | 2014 |
| d18O | BO14HTI02_d18O | BO14HTI01_d18O-BO14HTI02_d18O | 0.00 | Porites lutea | 28 | 1977 | 2005 | Pacific Ocean | Hon Tre Island, Vietnam | Bolton, A. | 2014 |
| d18O | LI99CLI01_d18O | LI99CLI01_d18O-WU14CLI01_d18O | 0.00 | Porites lobata | 24 | 1970 | 1994 | Pacific Ocean | Clipperton Island | Linsley, B. K. | 1999 |
| d18O | WU14CLI01_d18O | LI99CLI01_d18O-WU14CLI01_d18O | 0.00 | Porites lobata | 24 | 1970 | 1994 | Pacific Ocean | Clipperton Island | Wu, H. | 2014 |
| d18O | EV98KIR01_d18O | EV98KIR01_d18O-NU09KIR01_d18O | 18.47 | Porites sp. | 22 | 1972 | 1994 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | Evans, M. N. | 1998 |
| d18O | NU09KIR01_d18O | EV98KIR01_d18O-NU09KIR01_d18O | 18.47 | Porites sp. | 22 | 1972 | 1994 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | Nurhati, I. S. | 2009 |
| d18O | HE18COC01_d18O | HE18COC01_d18O-HE18COC02_d18O | 1.01 | Porites sp. | 22 | 1988 | 2010 | Indian Ocean | Cocos (Keeling) Islands, Australia | Hennekam, R. | 2018 |
| d18O | HE18COC02_d18O | HE18COC01_d18O-HE18COC02_d18O | 1.01 | Porites sp. | 22 | 1988 | 2010 | Indian Ocean | Cocos (Keeling) Islands, Australia | Hennekam, R. | 2018 |
| d18O | MC11KIR01_d18O | MC11KIR01_d18O-NU09KIR01_d18O | 18.47 | Porites sp. | 20 | 1978 | 1998 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | McGregor, H. V. | 2011 |
| d18O | NU09KIR01_d18O | MC11KIR01_d18O-NU09KIR01_d18O | 18.47 | Porites sp. | 20 | 1978 | 1998 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | Nurhati, I. S. | 2009 |
| d18O | OS13NGP01_d18O | OS13NGP01_d18O-OS13NLP01_d18O | 31.19 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Ngaragabel, Palau | Osborne, M. C. | 2013 |
| d18O | OS13NLP01_d18O | OS13NGP01_d18O-OS13NLP01_d18O | 31.19 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Ngeralang, Palau | Osborne, M. C. | 2013 |
| d18O | OS13NLP01_d18O | OS13NLP01_d18O-OS14RIP01_d18O | 47.16 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Ngeralang, Palau | Osborne, M. C. | 2013 |
| d18O | OS14RIP01_d18O | OS13NLP01_d18O-OS14RIP01_d18O | 47.16 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Rock Islands, Palau | Osborne, M. C. | 2014 |
| d18O | OS13NLP01_d18O | OS13NLP01_d18O-OS14UCP01_d18O | 53.77 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Ngeralang, Palau | Osborne, M. C. | 2013 |
| d18O | OS14UCP01_d18O | OS13NLP01_d18O-OS14UCP01_d18O | 53.77 | Porites lutea | 18 | 1990 | 2008 | Pacific Ocean | Ulong Channel, Palau | Osborne, M. C. | 2014 |
| d18O | EV98KIR01_d18O | EV98KIR01_d18O-MC11KIR01_d18O | 0.00 | Porites sp. | 16 | 1978 | 1994 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | Evans, M. N. | 1998 |
| d18O | MC11KIR01_d18O | EV98KIR01_d18O-MC11KIR01_d18O | 0.00 | Porites sp. | 16 | 1978 | 1994 | Pacific Ocean | Kiritimati (Christmas) Island, Republic of Kiribati | McGregor, H. V. | 2011 |
| d18O | AB20MEN05_d18O | AB20MEN05_d18O-AB20MEN06_d18O | 0.06 | Porites sp. | 15 | 1533 | 1548 | Indian Ocean | Mentawai Islands, Indonesia | Abram, N. | 2020 |
| d18O | AB20MEN06_d18O | AB20MEN05_d18O-AB20MEN06_d18O | 0.06 | Porites sp. | 15 | 1533 | 1548 | Indian Ocean | Mentawai Islands, Indonesia | Abram, N. | 2020 |
| d18O | NU11PAL01_d18O | NU11PAL01_d18O-SA19PAL01_d18O | 1.57 | Porites lutea | 14 | 1984 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Nurhati, I. S. | 2011 |
| d18O | SA19PAL01_d18O | NU11PAL01_d18O-SA19PAL01_d18O | 1.57 | Porites sp. | 14 | 1984 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| d18O | CO03PAL07_d18O | CO03PAL07_d18O-CO03PAL08_d18O | 0.00 | Porites lutea | 13 | 1654 | 1667 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | CO03PAL08_d18O | CO03PAL07_d18O-CO03PAL08_d18O | 0.00 | Porites lutea | 13 | 1654 | 1667 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Cobb, K. M. | 2003 |
| d18O | NU11PAL01_d18O | NU11PAL01_d18O-SA19PAL02_d18O | 1.57 | Porites lutea | 12 | 1986 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Nurhati, I. S. | 2011 |
| d18O | SA19PAL02_d18O | NU11PAL01_d18O-SA19PAL02_d18O | 1.57 | Porites sp. | 12 | 1986 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| d18O | SA19PAL01_d18O | SA19PAL01_d18O-SA19PAL02_d18O | 0.00 | Porites sp. | 12 | 1986 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| d18O | SA19PAL02_d18O | SA19PAL01_d18O-SA19PAL02_d18O | 0.00 | Porites sp. | 12 | 1986 | 1998 | Pacific Ocean | Palmyra Island, United States Minor Outlying Islands | Sayani, H. | 2019 |
| d18O | RE19GBR01_d18O | RE19GBR01_d18O-RE19GBR02_d18O | 26.34 | Porites sp. | 11 | 1979 | 1990 | Pacific Ocean | Eel Reef, Australia | Reed, E. V. | 2019 |
| d18O | RE19GBR02_d18O | RE19GBR01_d18O-RE19GBR02_d18O | 26.34 | Porites sp. | 11 | 1979 | 1990 | Pacific Ocean | Portland Roads, Australia | Reed, E. V. | 2019 |
| d18O | KR20SAR01_d18O | KR20SAR01_d18O-KR20SAR02_d18O | 9.82 | Porites sp. | 10 | 2006 | 2016 | Pacific Ocean | Sarawak, Malaysia | Krawczyk, H. | 2020 |
| d18O | KR20SAR02_d18O | KR20SAR01_d18O-KR20SAR02_d18O | 9.82 | Porites sp. | 10 | 2006 | 2016 | Pacific Ocean | Sarawak, Malaysia | Krawczyk, H. | 2020 |
