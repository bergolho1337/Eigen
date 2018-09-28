cmake_minimum_required (VERSION 3.0)
project (main)

set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(solver)

add_executable (main main.cpp)

add_dependencies(main solvers)

target_link_libraries (main solvers Eigen3::Eigen)