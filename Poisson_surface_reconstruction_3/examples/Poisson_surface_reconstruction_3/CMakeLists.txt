# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.15)
project(Poisson_surface_reconstruction_3_Examples)

# Find CGAL
find_package(CGAL REQUIRED)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  # Allow Windows 32bit applications to use up to 3GB of RAM
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

  # Print new compilation options
  message(
    STATUS
      "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(
    STATUS
      "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'"
  )
  message(
    STATUS
      "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'"
  )
  message(
    STATUS
      "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'"
  )
endif()

# Find Eigen3 (requires 3.1.0 or greater)
find_package(Eigen3 3.1.0)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  # Executables that require Eigen 3
  create_single_source_cgal_program("poisson_reconstruction_example.cpp")
  target_link_libraries(poisson_reconstruction_example
                        PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("poisson_reconstruction.cpp")
  target_link_libraries(poisson_reconstruction PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("poisson_reconstruction_function.cpp")
  target_link_libraries(poisson_reconstruction_function
                        PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("tutorial_example.cpp")
  target_link_libraries(tutorial_example PUBLIC CGAL::Eigen_support)
else()
  message(
    STATUS
      "NOTICE: The examples need Eigen 3.1 (or greater) will not be compiled.")
endif()
