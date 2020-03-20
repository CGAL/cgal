# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_deformation_Tests )



find_package(CGAL QUIET)

if ( CGAL_FOUND )

  find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
  if (EIGEN3_FOUND)
    create_single_source_cgal_program( "Cactus_deformation_session.cpp" )
    CGAL_target_use_Eigen(Cactus_deformation_session)
    create_single_source_cgal_program( "Cactus_performance_test.cpp" )
    CGAL_target_use_Eigen(Cactus_performance_test)
    create_single_source_cgal_program( "Symmetry_test.cpp" )
    CGAL_target_use_Eigen(Symmetry_test)

    find_package( OpenMesh QUIET )
    if ( OpenMesh_FOUND )
    include( UseOpenMesh )
      create_single_source_cgal_program( "Cactus_deformation_session_OpenMesh.cpp" )
      CGAL_target_use_Eigen(Cactus_deformation_session_OpenMesh)
      target_link_libraries( Cactus_deformation_session_OpenMesh PRIVATE ${OPENMESH_LIBRARIES} )
    else()
      message(STATUS "Example that use OpenMesh will not be compiled.")
    endif()
  else()
    message(STATUS "NOTICE: These tests require the Eigen library, version 3.2 or later and will not be compiled.")
  endif()
else()
  
    message(STATUS "NOTICE: These tests require the CGAL library, and will not be compiled.")
  
endif()

