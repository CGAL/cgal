# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( adding_handles.cpp )
  create_single_source_cgal_program( color.cpp )
  create_single_source_cgal_program( find_conflicts.cpp )
  create_single_source_cgal_program( hierarchy.cpp )
  create_single_source_cgal_program( linking_2d_and_3d.cpp )
  create_single_source_cgal_program( regular.cpp )
  create_single_source_cgal_program( simple.cpp )
  create_single_source_cgal_program( simplex.cpp )
  create_single_source_cgal_program( tds.cpp )

endif()
