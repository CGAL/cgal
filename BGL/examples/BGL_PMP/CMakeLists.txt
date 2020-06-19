project( BGL_PMP )


cmake_minimum_required(VERSION 2.8.11)

find_package(pmp)

if ( NOT pmp_FOUND )

  message(STATUS "This project requires the PMP library, and will not be compiled.")
  return()

endif()


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

find_package(Eigen3 3.2.0)
if (EIGEN3_FOUND)
  # Executables that require Eigen 3.2
  include( CGAL_Eigen_support )
endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()


create_single_source_cgal_program( "bgl_pmp.cpp" )

target_link_libraries(bgl_pmp PRIVATE pmp)
target_link_libraries(bgl_pmp PUBLIC CGAL::Eigen_support)
