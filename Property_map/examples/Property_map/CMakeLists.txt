cmake_minimum_required(VERSION 3.1...3.23)
project(Property_map_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("dynamic_properties.cpp")

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("custom_property_map.cpp")
  target_link_libraries(custom_property_map PUBLIC CGAL::Eigen3_support)
endif()
