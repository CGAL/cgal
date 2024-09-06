cmake_minimum_required(VERSION 3.1...3.15)

project( Hyperbolic_surface_triangulation_2_demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)


# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

if ( NOT CGAL_Qt6_FOUND  OR NOT Qt6_FOUND)

  message(STATUS "This project requires the Qt6 library, and will not be compiled.")
  return()

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# ui files, created with Qt Designer
qt6_wrap_ui(UIS drawing_window_description.ui)

add_executable( Hyperbolic_surface_triangulation_2_demo Hyperbolic_surface_triangulation_2_demo.cpp window.cpp ${UIS})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Hyperbolic_surface_triangulation_2_demo )

target_link_libraries(Hyperbolic_surface_triangulation_2_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets )

set(CMAKE_BUILD_TYPE "Release")

target_include_directories(Hyperbolic_surface_triangulation_2_demo PRIVATE ../../include/)
