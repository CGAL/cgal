# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_Examples )

cmake_minimum_required(VERSION 2.8.11)

# If you want to visualize a surface mesh, you can use the following viewer
# based on qt. Just uncomment the following line, plus the line find_package(CGALBasicViewer) below
# find_package(CGAL COMPONENTS Qt5)

# If you don't want to visualize, use the following line (otherwise comment it)
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  # find_package(CGALBasicViewer)

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "sm_derivation.cpp" )
  create_single_source_cgal_program( "sm_join.cpp" )
  if (USE_BASIC_VIEWER)
    target_link_libraries("sm_join" ${BASIC_VIEWER_LIBRARIES})
    qt5_use_modules("sm_join" ${BASIC_VIEWER_MODULES})
  endif(USE_BASIC_VIEWER)
  
  create_single_source_cgal_program( "sm_aabbtree.cpp" )
  create_single_source_cgal_program( "sm_bgl.cpp" )
  create_single_source_cgal_program( "sm_circulators.cpp" )
  create_single_source_cgal_program( "sm_do_intersect.cpp" )
  create_single_source_cgal_program( "sm_iterators.cpp" )
  create_single_source_cgal_program( "sm_kruskal.cpp" )
  create_single_source_cgal_program( "sm_memory.cpp" )
  create_single_source_cgal_program( "sm_properties.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

