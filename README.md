# omplib
Optical Model Potential LIBrary

This is primarily a library to be integrated into other projects, but can be used as a framework for standalone applications. For more information and references, build the documentation as shown below.

## quickstart

```
git clone git@github.com:beykyle/omplib.git
cd omplib
mkdir build
cd build 
cmake ..
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

## Python integration

This code will eventually include `Python` drivers using `numpy` and `xtensor`. 

For now, the only way to drive this code from `Python` is by adding the following `CMake` flag: `-DBUILD_SHARED_LIBS=True`. Then, the shared object file  `build/src/libomplib.so` will be created, which can interface with `Python` using the `ctypes` module. The `extern c` interface is provided below:
```
```

## building a standalone application
A (rather barebones) file with a main function lives in `exec/omplib.cpp`. To use `OMPLib` as a framework, copy the `exec` directory to a new directory (`my_exec`) of your choice within the main project directory. Change the name of `my_exec/omplib.cpp` (to e.g. `my_main.cpp`) as you see fit, and add the code for your application. Then replace `omplib.cpp` in `my_exec/CMakeLists.txt`, so it reads:

```
add_executable(my_exec_name my_main.cpp)
target_link_libraries(my_exec_name omplib)
```

Then, replace `add_subdirectory(exec)` with `add_subdirectory(my_exec)` in the top-level `CMakeLists.txt`. Finally, build as shown above.
