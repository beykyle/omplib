# cliutils.gitlab.io/modern-cmake/chapters/basics.html
# https://docs.travis-ci.com/user/languages/cpp/#cmake
cmake_minimum_required(VERSION 3.12.4)

# TODO this is a bit hacky
if(DEFINED PROJECT_NAME)
 set(SUBPROJECT ON)
else()
  set(SUBPROJECT OFF)
  project("OMPLib" VERSION 0.1 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
# valloric.github.io/YouCompleteMe/#c-family-semantic-completion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# stackoverflow.com/a/6595001/5101335
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif ()

find_package(Boost 1.7 REQUIRED)

find_package(Git REQUIRED)
find_package(Threads)

# stackoverflow.com/a/24414345
if (MSVC)
  add_compile_options(/permissive- /std:c++20)
  add_compile_definitions(_USE_MATH_DEFINES)
endif (MSVC)

# target_link_options() not available in CMake 3.12
if (!MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage")
endif (!MSVC)

### External projects
#####################
include(FetchContent)

FetchContent_Declare(
  json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
Message("fetching nlohmann/json")
FetchContent_MakeAvailable(json)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
Message("fetching Eigen")
FetchContent_MakeAvailable(Eigen)
set(EIGEN3_INCLUDE_DIR ${eigen_SOURCE_DIR})

## SUBDIRECTORIES
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(docs)
if (NOT SUBPROJECT)
  add_subdirectory(exec)
endif()
