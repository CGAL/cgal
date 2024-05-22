cmake_minimum_required(VERSION 3.1...3.15)

project( examples )

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

add_executable( example example.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS example)

target_link_libraries(example PRIVATE CGAL::CGAL )

set(CMAKE_BUILD_TYPE "Release")

#########################################################################
target_include_directories(example PRIVATE ../../include/)
