pybind11_add_module(omplibpy omplibpy.cpp)

target_link_libraries(omplibpy PRIVATE omplib)
target_include_directories(omplibpy PUBLIC
  ${CMAKE_SOURCE_DIR}/src ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

target_compile_definitions(omplibpy
  PRIVATE VERSION_INFO=${OMPLIBPY_VERSION_INFO})
