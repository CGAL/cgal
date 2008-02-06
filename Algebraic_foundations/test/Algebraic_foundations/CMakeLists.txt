# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( Algebraic_extension_traits.cpp )
  create_single_source_cgal_program( Algebraic_structure_traits.cpp )
  create_single_source_cgal_program( Coercion_traits.cpp )
  create_single_source_cgal_program( extended_euclidean_algorithm.cpp )
  create_single_source_cgal_program( Real_embeddable_traits.cpp )
  create_single_source_cgal_program( Scalar_factor_traits.cpp )

endif()
