# cliutils.gitlab.io/modern-cmake/chapters/basics.html
# https://docs.travis-ci.com/user/languages/cpp/#cmake
cmake_minimum_required(VERSION 3.12.4)
project("OMPLib" VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
# valloric.github.io/YouCompleteMe/#c-family-semantic-completion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# stackoverflow.com/a/6595001/5101335
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif ()

find_package (Git REQUIRED)
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
  json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
  DOWNLOAD_EXTRACT_TIMESTAMP true
  )
FetchContent_MakeAvailable(json)

## SUBDIRECTORIES
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(docs)
