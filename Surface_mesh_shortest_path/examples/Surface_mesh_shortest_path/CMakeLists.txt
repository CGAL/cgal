# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_shortest_path_ )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "shortest_path_sequence.cpp" )
  create_single_source_cgal_program( "shortest_paths_multiple_sources.cpp" )
  create_single_source_cgal_program( "shortest_paths_no_id.cpp" )
  create_single_source_cgal_program( "shortest_paths_with_id.cpp" )
  
  find_package( OpenMesh QUIET )
  if ( OpenMesh_FOUND )
    include( UseOpenMesh )
    create_single_source_cgal_program( "shortest_paths_OpenMesh.cpp" )
    target_link_libraries( shortest_paths_OpenMesh ${OPENMESH_LIBRARIES} )
  else()
    message(STATUS "Examples that use OpenMesh will not be compiled.")
  endif()
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

