add_executable(omp_run omplib.cpp)
target_link_libraries(omp_run PRIVATE nlohmann_json::nlohmann_json omplib)
