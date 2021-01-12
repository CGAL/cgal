# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.1...3.15)
project(Ridges_3_Examples)

find_package(CGAL REQUIRED)

# use either Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost QUIET COMPONENTS program_options)
  if(TARGET Boost::program_options)
    set(Boost_PROGRAM_OPTIONS_LIBRARY Boost::program_options)
  endif()
  if(Boost_PROGRAM_OPTIONS_FOUND)
    if(CGAL_AUTO_LINK_ENABLED)
      message(STATUS "Boost.ProgramOptions library: found")
    else()
      message(
        STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
    endif()
    add_definitions("-DCGAL_USE_BOOST_PROGRAM_OPTIONS")
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
  endif()

  create_single_source_cgal_program(Compute_Ridges_Umbilics.cpp)
  target_link_libraries(Compute_Ridges_Umbilics PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program(Ridges_Umbilics_SM.cpp)
  target_link_libraries(Ridges_Umbilics_SM PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program(Ridges_Umbilics_LCC.cpp)
  target_link_libraries(Ridges_Umbilics_LCC PUBLIC CGAL::Eigen_support)

else()

  message(
    STATUS
      "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled."
  )

endif()
