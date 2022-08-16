# cliutils.gitlab.io/modern-cmake/chapters/basics.html
# https://docs.travis-ci.com/user/languages/cpp/#cmake
cmake_minimum_required(VERSION 3.12.4)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
endif ()

find_package (Git)
if (GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe OUTPUT_VARIABLE output RESULT_VARIABLE result)
  string(COMPARE EQUAL ${result} "0" GIT_VERSION_FOUND)
  if (GIT_VERSION_FOUND)
    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" version "${output}")
  else(GIT_VERSION_FOUND)
    set(version "0.0.0")
    message("Git did not return a matching version number. You may have to `git fetch --tags` from remote.")
  endif(GIT_VERSION_FOUND)
endif (GIT_FOUND)

project("OMPLib" VERSION ${version} LANGUAGES CXX)

include(CTest)
find_package (Threads)

set(CMAKE_CXX_STANDARD 14)
# valloric.github.io/YouCompleteMe/#c-family-semantic-completion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# stackoverflow.com/a/6595001/5101335
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## copy over single header includes to build dir
# copy json
file(COPY json/single_include/nlohmann/json.hpp DESTINATION ../src/)

## SUBDIRECTORIES
add_subdirectory(src)
