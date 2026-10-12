# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.12...3.31)
project(Ridges_3_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    create_single_source_cgal_program(Compute_Ridges_Umbilics.cpp)
    target_link_libraries(Compute_Ridges_Umbilics PRIVATE CGAL::Eigen3_support)
    create_single_source_cgal_program(Ridges_Umbilics_SM.cpp)
    target_link_libraries(Ridges_Umbilics_SM PRIVATE CGAL::Eigen3_support)
    create_single_source_cgal_program(Ridges_Umbilics_LCC.cpp)
    target_link_libraries(Ridges_Umbilics_LCC PRIVATE CGAL::Eigen3_support)

    add_compile_definitions(CGAL_USE_BOOST_PROGRAM_OPTIONS)
    target_link_libraries(Compute_Ridges_Umbilics PRIVATE Boost::program_options)
    target_link_libraries(Ridges_Umbilics_SM PRIVATE Boost::program_options)
    target_link_libraries(Ridges_Umbilics_LCC PRIVATE Boost::program_options)

  else()
    message("NOTICE: This project requires Boost Program Options and will not be compiled.")
  endif()
else()
  message("NOTICE: This project requires Eigen 3.1 (or greater), and will not be compiled.")
endif()
