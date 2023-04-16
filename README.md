# ctmicro

> **Warning**
> The package is broken for Cantera 3.0, which is related to changes of 1D object instantiation. For guidance how to implement 1D objects, see [ctapp](https://github.com/ischoegl/ctapp).

The code defines a derived cantera class that implements heat transfer within a
channel with prescribed temperature (i.e. it allows for lateral heat losses,
making it a 1D implementation of non-adiabatic combustion)

**Citation:**

> Mohsen Ayoobi and Ingmar Schoegl, *Non-Catalytic Conversion of Glycerol to Syngas at Intermediate Temperatures: Numerical Methods with Detailed Chemistry*
> 2017, Fuel 195:190-200, [[DOI]](https://doi.org/10.1016/j.fuel.2017.01.065).


## Conda Environment

```
$ conda create -n ctmicro -c cantera/label/dev libcantera-devel cantera numpy scipy jupyter cython
$ conda activate ctmicro
```

## Usage

**Installation:** compile and install a *linked* version within the current python environment

```
$ python setup.py build_ext --inplace
$ pip install -e .
```

**Update:** run

```
$ git pull
$ python setup.py develop
```

**Unit Tests:** run

```
$ python setup.py test
```

**Uninstall:**

```
$ pip uninstall ctmicro
```
