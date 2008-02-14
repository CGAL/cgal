# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Surface_mesh_parameterization_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( Authalic_parameterization.cpp )
  create_single_source_cgal_program( Complete_parameterization_example.cpp )
  create_single_source_cgal_program( Mesh_cutting_parameterization.cpp )
  create_single_source_cgal_program( polyhedron_ex_parameterization.cpp )
  create_single_source_cgal_program( Simple_parameterization.cpp )
  create_single_source_cgal_program( Square_border_parameterization.cpp )
  create_single_source_cgal_program( Taucs_parameterization.cpp )

endif()
