# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Filtered_kernel_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
add_executable(bench_simple_comparisons bench_simple_comparisons.cpp)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  add_executable( bench_orientation_3 "orientation_3.cpp" )
  target_link_libraries(bench_orientation_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  add_executable( bench_comparisons "orientation_3.cpp" )
  target_link_libraries(bench_comparisons ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  set_property(
    TARGET bench_comparisons 
    APPEND
    PROPERTY COMPILE_DEFINITIONS ONLY_TEST_COMPARISONS )

  get_property(DEF TARGET bench_comparisons PROPERTY COMPILE_DEFINITIONS)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

