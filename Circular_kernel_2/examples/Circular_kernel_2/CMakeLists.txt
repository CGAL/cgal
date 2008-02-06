# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( arrangement_random_circles.cpp )
  create_single_source_cgal_program( arrangement_random_circles_segments.cpp )
  create_single_source_cgal_program( Circular_kernel_arrangement.cpp )
  create_single_source_cgal_program( Circular_kernel_arrangement_Line_arc.cpp )
  create_single_source_cgal_program( Circular_kernel_arrangement_Line_Circular.cpp )
  create_single_source_cgal_program( Exact_circular_kernel.cpp )
  create_single_source_cgal_program( Lazy_Circular_kernel_arrangement.cpp )
  create_single_source_cgal_program( Lazy_Circular_kernel_arrangement_Line_arc.cpp )
  create_single_source_cgal_program( Lazy_Circular_kernel_arrangement_Line_Circular.cpp )

endif()
