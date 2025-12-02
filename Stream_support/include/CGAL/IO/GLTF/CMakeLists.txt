# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
set(CMAKE_TOOLCHAIN_FILE "C:/Users/deman/vcpkg/scripts/buildsystems/vcpkg.cmake")
project( GLTF )

find_package(CGAL REQUIRED)

find_package(Eigen3 3.3.0 REQUIRED)
include(CGAL_Eigen3_support)

add_executable(test_read_GLTF test_read_GLTF.cpp tiny_gltf_impl.cpp)
target_link_libraries(test_read_GLTF PUBLIC CGAL::CGAL CGAL::Eigen3_support)