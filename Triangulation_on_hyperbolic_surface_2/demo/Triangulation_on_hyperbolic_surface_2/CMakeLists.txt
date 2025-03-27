cmake_minimum_required(VERSION 3.12...3.31)

project( Triangulation_on_hyperbolic_surface_2_Demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core Qt6)
find_package(Qt6 QUIET COMPONENTS Widgets)

if (NOT CGAL_Qt6_FOUND  OR NOT Qt6_FOUND)
  message(STATUS "This project requires the Qt6 library, and will not be compiled.")
  return()
endif()

# ui files, created with Qt Designer
qt6_wrap_ui(UIS drawing_window_description.ui)

add_executable( Triangulation_on_hyperbolic_surface_2_demo Triangulation_on_hyperbolic_surface_2_demo.cpp window.cpp ${UIS})
target_link_libraries(Triangulation_on_hyperbolic_surface_2_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets )
cgal_add_compilation_test(Triangulation_on_hyperbolic_surface_2_demo)
