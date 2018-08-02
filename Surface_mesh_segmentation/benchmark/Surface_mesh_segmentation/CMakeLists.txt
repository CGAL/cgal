# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_segmentation_Benchmarks )


cmake_minimum_required(VERSION 2.8.11)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# Activate concurrency option
option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)
if ( ACTIVATE_CONCURRENCY )
    find_package( TBB REQUIRED )

    if ( TBB_FOUND )
        include(${TBB_USE_FILE})
        list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    endif()
endif()

# include for local package
include_directories( BEFORE ../../include )



# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "approximate_convex_segmentation_benchmark.cpp" )
create_single_source_cgal_program( "segmentation_using_sdf_values_benchmark.cpp" )
