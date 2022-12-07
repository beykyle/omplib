if(NOT CMAKE_TESTING_ENABLED)
  Include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.0.1 # or a later release
  )
  Message("fetching Catch2")
  FetchContent_MakeAvailable(Catch2)


  # add new unit test source files here
  set(UNIT_TEST_SOURCES
    test_read_json.cpp
    test_channel.cpp
    test_yamaguchi.cpp
    test_gauss_legendre.cpp
    )
  # Add unit test input files here
  # Prefixes are stripped so duplicate filenames must not appear
  set(UNIT_TEST_INPUTS
    data/KDUQSamples/KD_default.json
    data/CHUQSamples/CH89_default.json
    data/WLHSamples/WLH_mean.json
    )

  list(TRANSFORM UNIT_TEST_INPUTS PREPEND ${CMAKE_SOURCE_DIR}/)
  add_custom_target(
    copy_ut_inputs ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${UNIT_TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}
    )

  # https://github.com/catchorg/Catch2/blob/devel/docs/slow-compiles.md
  add_executable(all_tests catch_main.cpp ${UNIT_TEST_SOURCES})
  target_link_libraries(all_tests PRIVATE
    nlohmann_json::nlohmann_json Catch2::Catch2WithMain omplib)

  add_custom_target(test
    COMMAND all_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Run Catch2 tests")
  add_dependencies(test all_tests copy_ut_inputs)

endif()


