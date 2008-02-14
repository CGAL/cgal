# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Intersections_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( tst04.cpp )
  create_single_source_cgal_program( tst05.cpp )
  create_single_source_cgal_program( tst08.cpp )
  create_single_source_cgal_program( tst09.cpp )
  create_single_source_cgal_program( tst10.cpp )
  create_single_source_cgal_program( tst11.cpp )
  create_single_source_cgal_program( tst50.cpp )
  create_single_source_cgal_program( tst51.cpp )
  create_single_source_cgal_program( tst52.cpp )
  create_single_source_cgal_program( tst53.cpp )
  create_single_source_cgal_program( tst54.cpp )
  create_single_source_cgal_program( tst56.cpp )
  create_single_source_cgal_program( tst57.cpp )
  create_single_source_cgal_program( tst58.cpp )
  create_single_source_cgal_program( tst59.cpp )

endif()
