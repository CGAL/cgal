# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Generator_example)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Boost COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  create_single_source_cgal_program("random_grid.cpp")
  create_single_source_cgal_program("random_disc_2.cpp")
  if(TARGET Boost::program_options)
    target_link_libraries(random_grid PRIVATE Boost::program_options)
    target_link_libraries(random_disc_2 PRIVATE Boost::program_options)
  else()
    target_link_libraries(random_grid PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
    target_link_libraries(random_disc_2 PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
  endif()
else()
  message("NOTICE: The benchmarks requires Boost Program Options, and will not be compiled.")
endif()
