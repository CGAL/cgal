# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_segmentation_Examples )


cmake_minimum_required(VERSION 3.1)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package( Boost )

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

# Activate concurrency option
option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)
if ( ACTIVATE_CONCURRENCY )
    find_package( TBB )

    if ( TBB_FOUND )
        include(${TBB_USE_FILE})
        list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    else()
        message(STATUS "TBB not found, examples of approximate convex decomposition will use the sequential version of the code.")
    endif()
endif()

# include for local package


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "sdf_values_example.cpp" )

create_single_source_cgal_program( "segmentation_from_sdf_values_example.cpp" )

create_single_source_cgal_program( "segmentation_via_sdf_values_example.cpp" )

create_single_source_cgal_program( "segmentation_with_facet_ids_example.cpp" )

create_single_source_cgal_program( "segmentation_from_sdf_values_SM_example.cpp")
create_single_source_cgal_program( "segmentation_from_sdf_values_LCC_example.cpp")

create_single_source_cgal_program( "extract_segmentation_into_mesh_example.cpp")

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "segmentation_from_sdf_values_OpenMesh_example.cpp" )
  target_link_libraries( segmentation_from_sdf_values_OpenMesh_example PRIVATE ${OPENMESH_LIBRARIES} )
endif()

create_single_source_cgal_program( "concavity_value_example.cpp" )
create_single_source_cgal_program( "concavity_value_SM_example.cpp" )

create_single_source_cgal_program( "approximate_convex_segmentation_example.cpp" )
create_single_source_cgal_program( "approximate_convex_segmentation_SM_example.cpp" )

create_single_source_cgal_program( "approximate_convex_segmentation_convex_hulls_example.cpp" )
