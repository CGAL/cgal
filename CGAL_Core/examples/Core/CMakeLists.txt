cmake_minimum_required(VERSION 3.1...3.15)
project(Core_Examples)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core)

if(NOT CGAL_Core_FOUND)

  message(
    STATUS
      "This project requires the CGAL_Core library, and will not be compiled.")
  return()

endif()

# Boost and its components
find_package(Boost)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package
include(${CGAL_USE_FILE})

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("delaunay.cpp")
