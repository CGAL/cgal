# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

set(Qt5_DIR /usr/local/Cellar/qt5/5.5.1_2/lib/cmake/Qt5)

project( Hyperbolic_triangulation_2 )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )

find_package(Qt5 COMPONENTS Core Gui Widgets)

include(${CGAL_USE_FILE})

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "Hyperbolic_Delaunay_triangulation_2_demo.cpp" )

qt5_use_modules(Hyperbolic_Delaunay_triangulation_2_demo Core Gui Widgets)

create_single_source_cgal_program( "Hyperbolic_Dirichlet_region_2_demo.cpp" )

qt5_use_modules(Hyperbolic_Dirichlet_region_2_demo Core Gui Widgets)

create_single_source_cgal_program( "Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp" )

qt5_use_modules(Periodic_2_Delaunay_hyperbolic_triangulation_2_demo Core Gui Widgets)

create_single_source_cgal_program( "Hyperbolic_translations_2_demo.cpp" )

qt5_use_modules(Hyperbolic_translations_2_demo Core Gui Widgets)




