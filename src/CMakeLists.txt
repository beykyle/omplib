file(GLOB_RECURSE SRC_CPP_FILES "*.cpp")
list(FILTER SRC_CPP_FILES EXCLUDE REGEX ".*main.cpp")

add_library(omplib_cpp ${SRC_CPP_FILES})
target_include_directories(omplib_cpp PUBLIC ${CMAKE_SOURCE_DIR}/src)

add_executable(omp main.cpp)
target_link_libraries(omp omplib_cpp ${CMAKE_THREAD_LIBS_INIT})
