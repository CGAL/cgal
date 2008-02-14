# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Voronoi_diagram_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( vda_ag.cpp )
  create_single_source_cgal_program( vda_dt.cpp )
  create_single_source_cgal_program( vda_pt.cpp )
  create_single_source_cgal_program( vda_rt.cpp )
  create_single_source_cgal_program( vda_sdg.cpp )

endif()
