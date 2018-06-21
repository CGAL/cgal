
project( Property_map_Tests )


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
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
  include( UseOpenMesh )
  add_definitions( -DCGAL_USE_OPENMESH )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "test_property_map.cpp" )

create_single_source_cgal_program( "dynamic_property_map.cpp" )

create_single_source_cgal_program( "dynamic_properties_test.cpp" )

create_single_source_cgal_program( "kernel_converter_properties_test.cpp" )

if(OpenMesh_FOUND)
        target_link_libraries( dynamic_properties_test PRIVATE ${OPENMESH_LIBRARIES} )
endif()

