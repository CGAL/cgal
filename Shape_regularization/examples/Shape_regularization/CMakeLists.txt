# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(Shape_regularization_Examples)

cmake_minimum_required(VERSION 3.1...3.15)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)
  message(STATUS "Found CGAL")

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Find OSQP library and headers.
  find_package(OSQP REQUIRED)
  if(OSQP_FOUND)
    message(STATUS "Found OSQP")

    set(osqp_targets
    regularize_15_segments
    regularize_100_segments_angles
    regularize_100_segments_offsets
    regularize_simple)

    foreach(osqp_target ${osqp_targets})
      create_single_source_cgal_program("${osqp_target}.cpp")
      if(TARGET ${osqp_target})
        target_link_libraries(${osqp_target} PUBLIC ${OSQP_LIBRARIES})
        target_compile_definitions(${osqp_target} PUBLIC -DCGAL_USE_OSQP)
      endif()
    endforeach()

    # Use Eigen.
    find_package(Eigen3 REQUIRED)
    include(CGAL_Eigen3_support)
    if(TARGET CGAL::Eigen3_support)
      message(STATUS "Found Eigen")
      create_single_source_cgal_program("regularize_real_data_2.cpp")
      target_link_libraries(regularize_real_data_2 PUBLIC CGAL::Eigen3_support ${OSQP_LIBRARIES})
      target_compile_definitions(regularize_real_data_2 PUBLIC -DCGAL_USE_OSQP)

    else()
      message(STATUS "NOTICE: Eigen was not found. Eigen examples won't be available.")
    endif()
  else()
    message(STATUS "NOTICE: OSQP was not found. OSQP examples won't be available.")
  endif()

  create_single_source_cgal_program("regularize_framework.cpp")
  create_single_source_cgal_program("regularize_planes.cpp")
  create_single_source_cgal_program("regularize_contour.cpp")
  create_single_source_cgal_program("regularize_closed_contour.cpp")
  create_single_source_cgal_program("regularize_open_contour.cpp")

else()
  message(WARNING "This program requires the CGAL library, and will not be compiled.")
endif()
