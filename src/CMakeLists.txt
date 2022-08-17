file(GLOB_RECURSE SRC_CPP_FILES "*.cpp")

add_library(omplib ${SRC_CPP_FILES})
target_include_directories(omplib PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_include_directories(omplib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(omplib PUBLIC
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate -fcoverage-mapping -O0 -Wall -Wextra -Wdocumentation>
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:Clang,GNU>>:-O3>)
target_link_options(omplib PUBLIC
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate>)

target_link_libraries(omplib PRIVATE nlohmann_json::nlohmann_json ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(omplib PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )

add_executable(omp omplib.cpp)
target_link_libraries(omp omplib)
