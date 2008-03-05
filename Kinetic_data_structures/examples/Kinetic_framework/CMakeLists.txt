# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kinetic_framework_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( adding_a_certificate.cpp )
  create_single_source_cgal_program( defining_a_simulation_traits.cpp )
  create_single_source_cgal_program( listener.cpp )
  create_single_source_cgal_program( pointer_queue.cpp )
  create_single_source_cgal_program( ref_counted.cpp )
  create_single_source_cgal_program( trivial_kds.cpp )

endif()
