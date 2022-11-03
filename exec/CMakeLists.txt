add_executable(omp_run omplib.cpp)
target_link_libraries(omp_run omplib nlohmann_json::nlohmann_json)
