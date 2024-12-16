# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Mesh_2_Demo)

find_package(CGAL REQUIRED)

# conform target
add_executable(conform conform.cpp)
target_link_libraries(conform CGAL::CGAL)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS conform)

# mesh target
add_executable(mesh mesh.cpp)
target_link_libraries(mesh CGAL::CGAL)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS mesh)

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
cgal_add_compilation_test(conform)
cgal_add_compilation_test(mesh)
