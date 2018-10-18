project( Property_map_Examples )


cmake_minimum_required(VERSION 3.1)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory


# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program( "dynamic_properties.cpp" )


