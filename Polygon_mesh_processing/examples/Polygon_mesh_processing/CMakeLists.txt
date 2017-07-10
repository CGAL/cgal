# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygon_mesh_processing_Examples )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

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

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

find_package( TBB )
if( TBB_FOUND )
  include( ${TBB_USE_FILE} )
  list( APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES} )
else()
  message( STATUS "NOTICE: Intel TBB was not found. hausdorff_distance_example will use sequential code." )
endif()
create_single_source_cgal_program( "hausdorff_distance_remeshing_example.cpp")
if (EIGEN3_FOUND)
  # Executables that require Eigen 3.2
  include( ${EIGEN3_USE_FILE} )

  create_single_source_cgal_program( "hole_filling_example.cpp" )
  create_single_source_cgal_program( "hole_filling_example_SM.cpp" )
  create_single_source_cgal_program( "refine_fair_example.cpp")
endif(EIGEN3_FOUND)

create_single_source_cgal_program( "self_intersections_example.cpp" )
create_single_source_cgal_program( "stitch_borders_example.cpp" )
create_single_source_cgal_program( "compute_normals_example_Polyhedron.cpp"  CXX_FEATURES cxx_range_for )
create_single_source_cgal_program( "compute_normals_example.cpp" CXX_FEATURES cxx_range_for cxx_auto_type )
create_single_source_cgal_program( "point_inside_example.cpp")
create_single_source_cgal_program( "triangulate_faces_example.cpp")
create_single_source_cgal_program( "connected_components_example.cpp")
create_single_source_cgal_program( "face_filtered_graph_example.cpp")
create_single_source_cgal_program( "polygon_soup_example.cpp")
create_single_source_cgal_program( "triangulate_polyline_example.cpp")
create_single_source_cgal_program( "mesh_slicer_example.cpp")
#create_single_source_cgal_program( "remove_degeneracies_example.cpp")
create_single_source_cgal_program( "isotropic_remeshing_example.cpp")
create_single_source_cgal_program( "surface_mesh_intersection.cpp")
create_single_source_cgal_program( "corefinement_SM.cpp")
create_single_source_cgal_program( "corefinement_consecutive_bool_op.cpp" )
create_single_source_cgal_program( "corefinement_difference_remeshed.cpp" )
create_single_source_cgal_program( "corefinement_mesh_union.cpp" )
create_single_source_cgal_program( "corefinement_polyhedron_union.cpp" )
create_single_source_cgal_program( "random_perturbation_SM_example.cpp" )
create_single_source_cgal_program( "corefinement_LCC.cpp")
create_single_source_cgal_program( "hole_filling_example_LCC.cpp" )

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



