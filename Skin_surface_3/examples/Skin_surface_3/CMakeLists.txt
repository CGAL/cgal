# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Skin_surface_3_example ) 

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL)




if (CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE include)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/NGHK_skin_surface_simple.cpp")
    # The two following .cpp files are in dont_submit, and thus not in the
    # internal releases.
    create_single_source_cgal_program( "NGHK_skin_surface_simple.cpp" )
    create_single_source_cgal_program( "NGHK_skin_surface_subdiv.cpp" )
  endif()
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
