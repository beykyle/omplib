file(GLOB_RECURSE TEST_CPP_FILES "*.cpp")

Include(FetchContent)
FetchContent_Declare(
  catch_download
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.3
  )
FetchContent_MakeAvailable(catch_download)

# add new unit test source files here
set(UNIT_TEST_SOURCES
  test_example.cpp
  )

# https://github.com/catchorg/Catch2/blob/devel/docs/slow-compiles.md
add_executable(run_tests catch2_main.cpp ${UNIT_TEST_SOURCES})
target_link_libraries(run_tests Catch2::Catch2 omplib)

# replace CTest with Catch2
add_custom_target(test
  COMMAND run_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Run Catch2 tests")
add_dependencies(test run_tests)
