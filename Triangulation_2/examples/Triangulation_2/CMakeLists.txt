# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( adding_handles.cpp )
  create_single_source_cgal_program( colored_face.cpp )
  create_single_source_cgal_program( constrained.cpp )
  create_single_source_cgal_program( constrained_hierarchy_plus.cpp )
  create_single_source_cgal_program( constrained_plus.cpp )
  create_single_source_cgal_program( hierarchy.cpp )
  create_single_source_cgal_program( regular.cpp )
  create_single_source_cgal_program( terrain.cpp )
  create_single_source_cgal_program( triangulation_prog1.cpp )
  create_single_source_cgal_program( voronoi.cpp )

endif()
