# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Visibility_complex_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( circle_bit.cpp )
  create_single_source_cgal_program( circle_vc.cpp )
  create_single_source_cgal_program( ellipse_bit.cpp )
  create_single_source_cgal_program( ellipse_vc.cpp )
  create_single_source_cgal_program( point_bit.cpp )
  create_single_source_cgal_program( point_sp.cpp )
  create_single_source_cgal_program( point_vc.cpp )
  create_single_source_cgal_program( poly_bit.cpp )
  create_single_source_cgal_program( poly_sp.cpp )
  create_single_source_cgal_program( poly_vc.cpp )
  create_single_source_cgal_program( segment_bit.cpp )
  create_single_source_cgal_program( segment_sp.cpp )
  create_single_source_cgal_program( segment_vc.cpp )

endif()
