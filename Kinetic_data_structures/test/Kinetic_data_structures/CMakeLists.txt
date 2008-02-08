# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( active_objects_tables.cpp )
  create_single_source_cgal_program( Delaunay_triangulation_2.cpp )
  create_single_source_cgal_program( Delaunay_triangulation_3.cpp )
  create_single_source_cgal_program( exact_kds.cpp )
  create_single_source_cgal_program( instantaneous_kernel.cpp )
  create_single_source_cgal_program( numbers.cpp )
  create_single_source_cgal_program( numeric_kds.cpp )
  create_single_source_cgal_program( random_kds.cpp )
  create_single_source_cgal_program( regular_triangulation_3.cpp )
  create_single_source_cgal_program( simulator.cpp )
  create_single_source_cgal_program( solvers.cpp )
  create_single_source_cgal_program( timings.cpp )

endif()
