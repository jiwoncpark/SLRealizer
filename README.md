# SLRealizer

Catalog-level simulation of LSST DM stack measurements of
gravitationally-lensed quasars.

In this package we are experimenting with predicting `Source`, `Object`,
`DIASource` and `DIAObject` properties directly, by realizing
model lenses as mixtures of Gaussians. Our goals are to

1. Be able to generate large-volume catalog level lens simulations, to enable catalog-level lens finder training;

2. Enable a model-based interpretation of the LSST catalogs, for use in lens finding and light curve extraction.


## Set-up and testing
For now, this is done using the DM stack-style setup from the `desc_package_template`: from `bash`
```
$ source <SLRealizer install directory>/setup/setup.sh
$ nosetests <SLRealizer install directory>
```
We'd like to switch to `setup.py` soon, as this is likely to be a pure python, stack-free analysis project.

## Demo

Watch this space!

## People
* [Phil Marshall](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@drphilmarshall) (SLAC)
* [Jenny Kim Marshall](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@jennykim1016 (Stanford)
* [Mike Baumer](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@mbaumer (SLAC)
* [Steve Kahn](https://github.com/LSSTDESC/SLRealizer/issues/new?body=@stevkahn (SLAC)


## License, etc.

This is research in progress, using open source software which is available for
re-use under the BSD license. If you make use of any of the ideas or code in
this repository in your own research, please do cite us as "(Marshall et al, in
prep.)" and provide a link to this repository - but before then, we hope you
will get in touch! Please do drop us a line via the hyperlinked contact names
above, or by [writing us an
issue](https://github.com/LSSTDESC/SLRealizer/issues/new).
