# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.1...3.20)
project(Ridges_3_Examples)

find_package(CGAL REQUIRED)

# use either Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost QUIET COMPONENTS program_options)

  if(Boost_PROGRAM_OPTIONS_FOUND)
    create_single_source_cgal_program(Compute_Ridges_Umbilics.cpp)
    target_link_libraries(Compute_Ridges_Umbilics PUBLIC CGAL::Eigen3_support)
    create_single_source_cgal_program(Ridges_Umbilics_SM.cpp)
    target_link_libraries(Ridges_Umbilics_SM PUBLIC CGAL::Eigen3_support)
    create_single_source_cgal_program(Ridges_Umbilics_LCC.cpp)
    target_link_libraries(Ridges_Umbilics_LCC PUBLIC CGAL::Eigen3_support)
    if(TARGET Boost::program_options)
      target_link_libraries(Compute_Ridges_Umbilics PRIVATE Boost::program_options)
      target_link_libraries(Ridges_Umbilics_SM PRIVATE Boost::program_options)
      target_link_libraries(Ridges_Umbilics_LCC PRIVATE Boost::program_options)
    else()
      target_link_libraries(Compute_Ridges_Umbilics PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
      target_link_libraries(Ridges_Umbilics_SM PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
      target_link_libraries(Ridges_Umbilics_LCC PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif()
  else()
    message(
      STATUS
        "NOTICE: This programs require Boost Program Options and will not be compiled."
    )
  endif()
else()
  message(
    STATUS
      "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
