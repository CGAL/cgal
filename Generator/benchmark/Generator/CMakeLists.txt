# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Generator_example)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

find_package(Boost REQUIRED program_options)
include_directories(${Boost_PROGRAM_OPTIONS_INCLUDE_DIR})
add_definitions("-DCGAL_USE_BOOST_PROGRAM_OPTIONS")
list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})

create_single_source_cgal_program("random_grid.cpp")
create_single_source_cgal_program("random_disc_2.cpp")
