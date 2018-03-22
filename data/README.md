# SLRealizer data

This folder contains a number of data files useful for experimenting with `SLRealizer`. You can also download the example "toy" outputs described in the concept DESC Note by [Kim et al (2018)](https://github.com/LSSTDESC/SLRealizer/tree/issue/17/desc-note/doc/slrealizer-concept-note/main.pdf).

## Observation History

* Twinkles visit list (919 rows, 40kb): [twinkles_observation_history.csv](twinkles_observation_history.csv)

## Example outputs

These are large files, stored on dropbox. You can obtain them as follows:

* Kim et al (2018) Toy Object Catalog (2234 rows, 1.6Mb):
```
wget -O toy_object_catalog.csv "https://www.dropbox.com/s/ob51rxjexzuervl/toy_object_catalog.csv?dl=0"
```
* Kim et al (2018) Toy Source Catalog (473596 rows, 87Mb):
```
wget -O toy_source_catalog.csv "https://www.dropbox.com/s/muqui8eu3kxox2l/toy_source_catalog.csv?dl=0"
```

## Non-Lenses for Machine Learning Tests

* SDSS bright object catalog (96594 rows, 59Mb, see [Kim et al (2018)](https://github.com/LSSTDESC/SLRealizer/tree/issue/17/desc-note/doc/slrealizer-concept-note/main.pdf) for details): [sdss_object.csv](sdss_object.csv)
