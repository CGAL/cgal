# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_2_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

set(CMAKE_AUTOMOC TRUE)
set(CMAKE_AUTOUIC TRUE)
set(CMAKE_AUTORCC TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

find_package(CGAL COMPONENTS Qt5)
find_package(Qt5 QUIET COMPONENTS Widgets)

if ( NOT CGAL_FOUND OR NOT CGAL_Qt5_FOUND OR NOT Qt5_FOUND )
  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")
  return()
endif()

add_definitions(-DQT_NO_KEYWORDS)

#--------------------------------
# The "constrained Delaunay" demo: Constrained_Delaunay_triangulation_2
#--------------------------------

# The executable itself.
add_executable( Constrained_Delaunay_triangulation_2
  Constrained_Delaunay_triangulation_2.cpp)
target_link_libraries( Constrained_Delaunay_triangulation_2
  PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)
target_include_directories( Constrained_Delaunay_triangulation_2
  PRIVATE ./include)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Constrained_Delaunay_triangulation_2)

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------

# The executable itself.
add_executable  ( Delaunay_triangulation_2
  Delaunay_triangulation_2.cpp)
target_link_libraries( Delaunay_triangulation_2
  PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_2 )

#--------------------------------
# The "Regular" demo: Regular_triangulation_2
#--------------------------------

# The executable itself.
add_executable  ( Regular_triangulation_2
  Regular_triangulation_2.cpp)
target_link_libraries(Regular_triangulation_2
  PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Regular_triangulation_2 )

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
foreach(target
    Constrained_Delaunay_triangulation_2
    Delaunay_triangulation_2
    Regular_triangulation_2)
  cgal_add_compilation_test(${target})
endforeach()

