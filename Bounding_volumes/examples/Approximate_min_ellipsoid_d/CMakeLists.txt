cmake_minimum_required(VERSION 3.1...3.15)
project( Approximate_min_ellipsoid_d_Examples )


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # Use Eigen
  find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
  
  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    if(NOT (${cppfile} STREQUAL "ellipsoid.cpp") OR EIGEN3_FOUND)
      create_single_source_cgal_program( "${cppfile}" )
      if (EIGEN3_FOUND)
        get_filename_component(target ${cppfile} NAME_WE)
        CGAL_target_use_Eigen(${target})
      endif()
    endif()
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

