# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Skin_surface_3_example ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL)




if (CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE include)

# With iterator debugging of VC++this example takes forever
  add_definitions( -D_HAS_ITERATOR_DEBUGGING=0 )

  create_single_source_cgal_program( "NGHK_skin_surface_simple.cpp" )
  create_single_source_cgal_program( "NGHK_skin_surface_subdiv.cpp" )
  create_single_source_cgal_program( "skin_surface_simple.cpp" )
  create_single_source_cgal_program( "skin_surface_subdiv.cpp" )
  create_single_source_cgal_program( "skin_surface_subdiv_with_normals.cpp" )
  create_single_source_cgal_program( "union_of_balls_simple.cpp" )
  create_single_source_cgal_program( "union_of_balls_subdiv.cpp" )  
  
  find_package(ESBTL)
  if(ESBTL_FOUND)
    include_directories( ${ESBTL_INCLUDE_DIR} )
    create_single_source_cgal_program( "skin_surface_pdb_reader.cpp" )    
  else(ESBTL_FOUND)
    message(STATUS "NOTICE: skin_surface_pdb_reader.cpp requires ESBTL library, and will not be compiled.")  
  endif(ESBTL_FOUND)

else(CGAL_FOUND)

  message(STATUS "NOTICE: This directory requires CGAL library, and will not be compiled.")

endif(CGAL_FOUND)
