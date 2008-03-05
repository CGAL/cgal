# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Circular_kernel_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_Circular_kernel.cpp )
  create_single_source_cgal_program( test_Exact_circular_kernel.cpp )
  create_single_source_cgal_program( test_Filtered_bbox_circular_kernel.cpp )
  create_single_source_cgal_program( test_Line_arc.cpp )

endif()
