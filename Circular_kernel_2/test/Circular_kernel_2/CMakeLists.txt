# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Circular_kernel_2_Tests)

find_package(CGAL REQUIRED)

include_directories(BEFORE include)

# The following `include_directories` is used in the git layout
include_directories(BEFORE ../../../Kernel_23/test/Kernel_23/include)

# The following `include_directories` is used in the layout of the
# internal releases tarballs
include_directories(BEFORE ../Kernel_23/include)

create_single_source_cgal_program("test_Circular_kernel.cpp")
create_single_source_cgal_program("test_Exact_circular_kernel.cpp")
create_single_source_cgal_program("test_Filtered_bbox_circular_kernel.cpp")
create_single_source_cgal_program("test_Line_arc.cpp")
create_single_source_cgal_program("test_Circular_kernel_basic.cpp")
create_single_source_cgal_program("test_Exact_circular_kernel_basic.cpp")

if(CGAL_ENABLE_TESTING)
  set_tests_properties(
    "execution   of  test_Circular_kernel_basic"
    "execution   of  test_Exact_circular_kernel_basic"
    PROPERTIES RESOURCE_LOCK Circular_kernel_basic)
endif()
