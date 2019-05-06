# ctmicro

The code defines a derived cantera class that implements heat transfer within a
channel with prescribed temperature (i.e. it allows for lateral heat losses,
making it a 1D implementation of non-adiabatic combustion)

**Citation:**

> Mohsen Ayoobi and Ingmar Schoegl, *Non-Catalytic Conversion of Glycerol to Syngas at Intermediate Temperatures: Numerical Methods with Detailed Chemistry*
> 2017, Fuel 195:190-200, [DOI](https://doi.org/10.1016/j.fuel.2017.01.065)


## Conda Environment

```
$ conda create -n ctmicro -c cantera/label/dev libcantera-devel cantera numpy scipy jupyter cython
$ conda activate ctmicro
```

## Usage

Within your environment, compile and install package

```
$ python setup.py build_ext --inplace
$ pip install -e .
```
