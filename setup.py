from setuptools import setup, find_packages

setup(
    name='slrealizer',
    version='0.2.1',
    author='Phil Marshall',
    author_email='dr.phil.marshall@gmail.com',
    package_dir={'':'slrealizer'},
    packages=find_packages('slrealizer'),
    license='LICENSE.md',
    description='Catalog-level realization of simulated gravitational lenses.',
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url='https://github.com/LSSTDESC/SLRealizer',
    install_requires=[
        "pip==9.0.3",
        "numpy>=1.13",
        "future>=0.15",
        "matplotlib>=2.2.2",
        "astropy=2.0",
        "pandas>=0.20",
        "om10"],
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python'],
    keywords='physics'
)
