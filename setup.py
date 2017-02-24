from setuptools import setup

setup(# package information
      name="slrealizer",
      version="0.0",
      author="Phil Marshall",
      author_email="dr.phil.marshall@gmail.com",
      description="Catalog-level realization of simulated gravitational lenses",
      long_description=open("README.md").read(),
      url="https://github.com/LSSTDESC/SLRealizer",
      packages=['desc.slrealizer'],
      package_dir={'desc.slrealizer':'python/desc'},
      include_package_data=True,
      package_data={'slrealizer': ['data/*']},
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: BSD License",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
      ],
      install_requires=["numpy", "scipy", "matplotlib", "astropy", "om10"],
)

# The above code does not work. `python setup.py develop` completes with
# no complaints, but the package does not import. Here's the log:

# $ python setup.py develop
# running develop
# running egg_info
# writing requirements to slrealizer.egg-info/requires.txt
# writing slrealizer.egg-info/PKG-INFO
# writing top-level names to slrealizer.egg-info/top_level.txt
# writing dependency_links to slrealizer.egg-info/dependency_links.txt
# writing manifest file 'slrealizer.egg-info/SOURCES.txt'
# running build_ext
# Creating /Users/pjm/miniconda2/envs/lsst/lib/python2.7/site-packages/slrealizer.egg-link (link to .)
# Adding slrealizer 0.0 to easy-install.pth file
# ...
# Using /Users/pjm/miniconda2/envs/lsst/lib/python2.7/site-packages
# Finished processing dependencies for slrealizer==0.0
#
# $ python
# Python 2.7.12 |Continuum Analytics, Inc.| (default, Jul  2 2016, 17:43:17)
# [GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)] on darwin
# Type "help", "copyright", "credits" or "license" for more information.
# Anaconda is brought to you by Continuum Analytics.
# Please check out: http://continuum.io/thanks and https://anaconda.org
# >>> import desc.slrealizer
# Traceback (most recent call last):
#   File "<stdin>", line 1, in <module>
# ImportError: No module named slrealizer
# >>> import slrealizer
# Traceback (most recent call last):
#   File "<stdin>", line 1, in <module>
# ImportError: No module named slrealizer
# >>> import desc
# >>> dir(desc)
# ['__builtins__', '__doc__', '__file__', '__name__', '__package__', '__path__', 'pkgutil']
# >>> desc.__package__
# 'desc'

# So 'desc' is setup as a package, but desc.slrealizer is not.
