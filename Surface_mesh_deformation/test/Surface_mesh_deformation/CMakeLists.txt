# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_deformation_test )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")
  find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "Cactus_deformation_session.cpp" )
    create_single_source_cgal_program( "Cactus_performance_test.cpp" )
    create_single_source_cgal_program( "Symmetry_test.cpp" )

    find_package( OpenMesh QUIET )
    if ( OpenMesh_FOUND )
    include( UseOpenMesh )
      create_single_source_cgal_program( "Cactus_deformation_session_OpenMesh.cpp" )
      target_link_libraries( Cactus_deformation_session_OpenMesh ${OPENMESH_LIBRARIES} )
    else()
      message(STATUS "Example that use OpenMesh will not be compiled.")
    endif()
  else()
    message(STATUS "NOTICE: These tests require the Eigen library, version 3.2 or later and will not be compiled.")
  endif()
else()
  
    message(STATUS "NOTICE: These tests require the CGAL library, and will not be compiled.")
  
endif()

