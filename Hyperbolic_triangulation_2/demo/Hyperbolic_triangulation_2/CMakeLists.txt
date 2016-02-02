# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Hyperbolic_triangulation_2 )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5Core Qt5Giu Qt5Widgets )

# Find the QtWidgets library
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5Widgets)

# Add the include directories for the Qt 5 Widgets module to
# the compile lines.
include_directories(${Qt5Core_INCLUDE_DIRS})
include_directories(${Qt5Gui_INCLUDE_DIRS})
include_directories(${Qt5Widgets_INCLUDE_DIRS})

# Use the compile definitions defined in the Qt 5 Widgets module
add_definitions(${Qt5Core_DEFINITIONS})
add_definitions(${Qt5Gui_DEFINITIONS})
add_definitions(${Qt5Widgets_DEFINITIONS})

# Add compiler flags for building executables (-fPIE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Gui_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

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


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )


qt5_generate_moc( Hyperbolic_Delaunay_triangulation_2_demo.cpp Hyperbolic_Delaunay_triangulation_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_Delaunay_triangulation_2_demo.cpp" )
target_link_libraries(Hyperbolic_Delaunay_triangulation_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################


qt5_generate_moc( Hyperbolic_Dirichlet_region_2_demo.cpp Hyperbolic_Dirichlet_region_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_Dirichlet_region_2_demo.cpp" )
target_link_libraries(Hyperbolic_Dirichlet_region_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################


qt5_generate_moc( Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.moc )
create_single_source_cgal_program( "Periodic_2_Delaunay_hyperbolic_triangulation_2_demo.cpp" )
target_link_libraries(Periodic_2_Delaunay_hyperbolic_triangulation_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)

################

qt5_generate_moc( Hyperbolic_translations_2_demo.cpp Hyperbolic_translations_2_demo.moc )
create_single_source_cgal_program( "Hyperbolic_translations_2_demo.cpp" )
target_link_libraries(Hyperbolic_translations_2_demo CGAL CGAL_Qt5 Qt5::Core Qt5::Gui Qt5::Widgets)














