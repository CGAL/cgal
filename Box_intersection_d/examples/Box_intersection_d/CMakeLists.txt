# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Box_intersection_d_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( box_grid.cpp )
  create_single_source_cgal_program( custom_box_grid.cpp )
  create_single_source_cgal_program( minimal.cpp )
  create_single_source_cgal_program( minimal_self.cpp )
  create_single_source_cgal_program( proximity_custom_box_traits.cpp )
  create_single_source_cgal_program( triangle_self_intersect.cpp )
  create_single_source_cgal_program( triangle_self_intersect_pointers.cpp )

endif()
