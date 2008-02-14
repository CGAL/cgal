# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Algebraic_kernel_d_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( Algebraic_curve_kernel_2.cpp )
  create_single_source_cgal_program( Algebraic_kernel_1.cpp )
  create_single_source_cgal_program( Algebraic_real_pure.cpp )
  create_single_source_cgal_program( Bitstream_descartes.cpp )
  create_single_source_cgal_program( Descartes.cpp )

endif()
