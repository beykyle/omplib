# omplib
Optical Model Potential LIBrary

This is a library of microscopic and phenomenological potentials relevant for nuclear reactions, with a primary focus on nucleon-nuclear scattering using optical potentials. It is equipped with an R-Matrix solver on a Lagrange-Legendre mesh for calculating scattering matrix elements and cross sections. Applications include fitting and uncertainty quantification of optical model parameters, and Hauser-Feshbach inelastic scattering calculations. 

The goals of this project are to provide a fast, flexible and high-fidelity solver for scattering problems, including non-local potentials, coupled-channels, and more; that is easy to plug into existing physics packages. This is primarily built as a library to be integrated into other projects, but can be used as a framework for standalone applications. For more information and extensive references to the scientific literature, build the documentation as shown below.

## quickstart

```
git clone git@github.com:beykyle/omplib.git
cd omplib
mkdir build
cd build 
cmake -DCMAKE_BUILD_TYPE=Release .. 
make 
make test
make docs
```

## integration

OMPLib supports CMake integration using [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) with static linkage. Add the following to your `CMakeLists.txt`:

```
FetchContent_Declare(
  omplib GIT_REPOSITORY git@github.com:beykyle/omplib.git
  GIT_TAG "origin/main"
  )
message("Fetching beykyle/omplib")
FetchContent_MakeAvailable(omplib)
```

Now you can `#include` files like `"potential/params.hpp"` into your project, as long as you make `omplib` a dependency of the relevant target.

## python integration

`OMPLib` uses `pybind11` with `xtensor-python` to provide a `Python` module, `omplibpy`, with functionality to build potentials and drive the solvers in `OMPLib`. See the list of [dependencies](https://github.com/beykyle/omplib/blob/main/README.md#install-yourself-but-only-if-you-want-to-use-the-python-module-omplibpy) below for building the module. To manually build the `cpython` shared library for `OMPLib`, compile with the flag `-DBUILD_PY_MODULE=On`, with everything else otherwise unchanged as in [quickstart](https://github.com/beykyle/omplib/blob/main/README.md#quickstart). However, the best way to build the code and install the module all at once is using the `setup.py` script, e.g., buy running:

```
python setup.py install
```

or, 

```
pip install .
```

in the main project directory. As long as all the aforementioned dependencies exist in your `Python` environment, `omplibpy` will be built and installed, so you can `import` it, like so:
```
import omplibpy as omp

# 0.1 Mev neutron incident on 48-Ca, using WLH microscopic optical potential
erg_cms_Mev = 0.1
A = 20
Z = 48 
xs_tot, xs_rxn = omp.wlh_xs_n(0.1, A, Z)
```
See `examples/` for more simple `Python` scripts using the module.

## dependencies

### handled by `CMake` (just run `cmake ..`, kick your feet up, and relax):
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [nlohmann/json](https://github.com/nlohmann/json)

### install yourself (sorry):
- [Boost](https://www.boost.org/) 1.7+

### install yourself, but only if you want to use the `Python` module, `omplibpy`:
- [Python](https://www.python.org/) 3.7+
- [numpy](https://numpy.org/)
- [xtensor](https://xtensor.readthedocs.io/en/latest/)
- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)
- [xtensor-python](https://xtensor-python.readthedocs.io/en/latest/index.html)

It is highly recomended to use use a package, dependency and environment manager like [mamba](https://mamba.readthedocs.io/en/latest/) or [conda](https://docs.conda.io/en/latest/), especially if you want to use `omplibpy`. Then, setting up an environment to run `OMPLib` from `Python` is as easy as (e.g. using `mamba`):

```
mamba create -n omp boost numpy pybind11 xtensor xtensor-python
mamba activate omp
```

Now, within the environment `omp`, you can `pip install` as above to build `OMPLib` and install the `omplibpy` module.

## configuration

Currently, the primary application for `OMPLib` is uncertainty quantification for simple spherical potentials. For that reason, for single-channel calculations, it is set up to solve for the R-Matrix without heap allocations, using `Eigen` statically sized matrices, for speed. The file `src/util/config.hpp` holds important compile time constants, such as `NBASIS`; the number of basis functions to use in the R-Matrix solver. The resulting matrix for a single channel calculation will be represented as an `Eigen::Matrix<std::complex<double>,NBASIS,NBASIS>`. Inversion of these matrices is expected to dominate runtime for typical problems.
