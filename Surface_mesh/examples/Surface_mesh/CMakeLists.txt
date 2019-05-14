# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_Examples )

cmake_minimum_required(VERSION 2.8.11)

if(POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  if(NOT POLICY CMP0070)
    cmake_policy(SET CMP0053 OLD)
  else()
    cmake_policy(SET CMP0053 NEW)
  endif()
endif()
if(POLICY CMP0071)
    cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "sm_derivation.cpp" )
  create_single_source_cgal_program( "sm_join.cpp" )
  create_single_source_cgal_program( "sm_aabbtree.cpp" )
  create_single_source_cgal_program( "sm_bgl.cpp" )
  create_single_source_cgal_program( "sm_circulators.cpp" )
  create_single_source_cgal_program( "sm_do_intersect.cpp" )
  create_single_source_cgal_program( "sm_iterators.cpp" )
  create_single_source_cgal_program( "sm_kruskal.cpp" )
  create_single_source_cgal_program( "sm_memory.cpp" )
  create_single_source_cgal_program( "sm_properties.cpp" )

  create_single_source_cgal_program("draw_surface_mesh.cpp")
  if(CGAL_Qt5_FOUND )
    target_link_libraries(draw_surface_mesh PUBLIC CGAL::CGAL_Qt5)
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

