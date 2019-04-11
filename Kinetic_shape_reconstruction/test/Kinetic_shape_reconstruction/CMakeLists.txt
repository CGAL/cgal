# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( test_cmake )

cmake_minimum_required(VERSION 2.8.11)

###############################################################################
# TARGETS
###############################################################################
set(targets
  kinetic_2d_stress_test
)


###############################################################################
# GENERATED PART
###############################################################################

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

include( CGAL_CreateSingleSourceCGALProgram )

# Libraries and flags
set(project_linked_libraries)
set(project_compilation_definitions)

# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Creating targets with correct libraries and flags
foreach(target ${targets})
  create_single_source_cgal_program( "${target}.cpp" )
  if(TARGET ${target})
    target_link_libraries(${target} PUBLIC ${project_linked_libraries})
    target_compile_definitions(${target} PUBLIC ${project_compilation_definitions})
  endif()
endforeach()

