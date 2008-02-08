# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( tst321.cpp )
  create_single_source_cgal_program( tst322.cpp )
  create_single_source_cgal_program( tst323.cpp )
  create_single_source_cgal_program( tst324.cpp )
  create_single_source_cgal_program( tst325.cpp )
  create_single_source_cgal_program( tst326.cpp )
  create_single_source_cgal_program( tst327.cpp )
  create_single_source_cgal_program( tst328.cpp )
  create_single_source_cgal_program( tst329.cpp )
  create_single_source_cgal_program( tst330.cpp )
  create_single_source_cgal_program( tst331.cpp )
  create_single_source_cgal_program( tst332.cpp )
  create_single_source_cgal_program( tst333.cpp )

endif()
