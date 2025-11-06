# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.31)
project(Poisson_surface_reconstruction_3_Tests)

# Find CGAL
find_package(CGAL REQUIRED)
find_package(TBB QUIET)
include(CGAL_TBB_support)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  # Allow Windows 32bit applications to use up to 3GB of RAM
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

  # Prints new compilation options
  message(STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'")
  message(STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'")
  message(STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'")
endif()

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  # Executables that require Eigen 3.1
  if(TARGET CGAL::TBB_support)
    create_single_source_cgal_program("poisson_reconstruction_test_surface_mesher.cpp")
    target_link_libraries(poisson_reconstruction_test_surface_mesher PRIVATE CGAL::Eigen3_support CGAL::TBB_support)

    create_single_source_cgal_program("poisson_reconstruction_test_mesh_3.cpp")
    target_link_libraries(poisson_reconstruction_test_mesh_3 PRIVATE CGAL::Eigen3_support CGAL::TBB_support)

    create_single_source_cgal_program("compare_mesh_3_vs_Poisson_implicit_surface_3.cpp")
    target_link_libraries(compare_mesh_3_vs_Poisson_implicit_surface_3  PRIVATE CGAL::Eigen3_support CGAL::TBB_support)
  else()
    create_single_source_cgal_program("poisson_reconstruction_test_surface_mesher.cpp")
    target_link_libraries(poisson_reconstruction_test_surface_mesher PRIVATE CGAL::Eigen3_support)

    create_single_source_cgal_program("poisson_reconstruction_test_mesh_3.cpp")
    target_link_libraries(poisson_reconstruction_test_mesh_3 PRIVATE CGAL::Eigen3_support)

    create_single_source_cgal_program("compare_mesh_3_vs_Poisson_implicit_surface_3.cpp")
    target_link_libraries(compare_mesh_3_vs_Poisson_implicit_surface_3  PRIVATE CGAL::Eigen3_support)
  endif()
else()
  message("NOTICE: Tests in this directory require Eigen 3.1 (or greater), and will not be compiled.")
endif()

