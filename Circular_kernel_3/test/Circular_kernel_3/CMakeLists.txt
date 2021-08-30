# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Circular_kernel_3_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

include_directories(BEFORE include)

# The following `include_directories` is used in the git layout
include_directories(BEFORE ../../../Kernel_23/test/Kernel_23/include)

# The following `include_directories` is used in the layout of the
# internal releases tarballs
include_directories(BEFORE ../Kernel_23/include)

create_single_source_cgal_program("test_Spherical_kernel.cpp")
create_single_source_cgal_program("test_Spherical_kernel_basics.cpp")
create_single_source_cgal_program("test_Exact_spherical_kernel.cpp")
create_single_source_cgal_program("test_Lazy_Spherical_kernel.cpp")
create_single_source_cgal_program("test_Lazy_spherical_kernel_basics.cpp")
create_single_source_cgal_program("test_Spherical_kernel_with_core.cpp")
