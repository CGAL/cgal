# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygon_mesh_processing )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
 
  if (EIGEN3_FOUND)
    # Executables that require Eigen 3.1
    include( ${EIGEN3_USE_FILE} )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "self_intersections_example.cpp" )
create_single_source_cgal_program( "stitch_borders_example.cpp" )
create_single_source_cgal_program( "compute_normals_example.cpp" )
create_single_source_cgal_program( "point_inside_example.cpp")
create_single_source_cgal_program( "triangulate_faces_example.cpp")
create_single_source_cgal_program( "connected_components_example.cpp")
create_single_source_cgal_program( "polygon_soup_example.cpp")
create_single_source_cgal_program( "triangulate_polyline_example.cpp")
create_single_source_cgal_program( "mesh_slicer_example.cpp")
#create_single_source_cgal_program( "remove_degeneracies_example.cpp")

if(NOT (${EIGEN3_VERSION} VERSION_LESS 3.2.0))
create_single_source_cgal_program( "hole_filling_example.cpp" )
create_single_source_cgal_program( "hole_filling_example_SM.cpp" )
create_single_source_cgal_program( "refine_fair_example.cpp")
else()
  message(STATUS "NOTICE: Some examples require Eigen 3.2 (or greater) and will not be compiled.")  
endif()

if(OpenMesh_FOUND)
create_single_source_cgal_program( "compute_normals_example_OM.cpp" )
target_link_libraries( compute_normals_example_OM ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "hole_filling_example_OM.cpp" )
target_link_libraries( hole_filling_example_OM ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "point_inside_example_OM.cpp")
target_link_libraries( point_inside_example_OM ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "stitch_borders_example_OM.cpp" )
target_link_libraries( stitch_borders_example_OM ${OPENMESH_LIBRARIES} )

#create_single_source_cgal_program( "remove_degeneracies_example_OM.cpp")
#target_link_libraries( remove_degeneracies_example_OM ${OPENMESH_LIBRARIES} )

create_single_source_cgal_program( "triangulate_faces_example_OM.cpp")
target_link_libraries( triangulate_faces_example_OM ${OPENMESH_LIBRARIES} )
endif(OpenMesh_FOUND)

else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Some examples require Eigen 3.1 (or greater) and will not be compiled.")  
endif(EIGEN3_FOUND)

