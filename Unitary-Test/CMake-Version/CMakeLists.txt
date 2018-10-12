cmake_minimum_required (VERSION 3.0)
project (solver)

set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(solvers)

add_executable (solver solver.cpp)

add_dependencies(solver solvers)

target_link_libraries (solver solvers Eigen3::Eigen)