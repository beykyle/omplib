file(GLOB_RECURSE SRC_CPP_FILES "*.cpp")

add_library(omplib ${SRC_CPP_FILES})
target_include_directories(omplib PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_include_directories(omplib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

target_compile_options(omplib PUBLIC
  $<$<AND:$<CONFIG:Debug>>:
    -O0 -Wall -Wextra -g>
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:Clang,GNU>>:
   -O3>
)

target_link_libraries(omplib PUBLIC
  Eigen3::Eigen nlohmann_json::nlohmann_json ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(omplib PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
