cmake_minimum_required(VERSION 3.1...3.15)

project( demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)


# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core Qt5)
include(${CGAL_USE_FILE})
find_package(Qt5 QUIET COMPONENTS Widgets)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

if ( NOT CGAL_Qt5_FOUND )

  message(STATUS "This project requires the Qt5 library, and will not be compiled.")
  return()

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# ui files, created with Qt Designer
qt5_wrap_ui(UIS drawing_window_description.ui)

add_executable( demo demo.cpp window.cpp ${UIS})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS demo )

target_link_libraries(demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets )

set(CMAKE_BUILD_TYPE "Release")

target_include_directories(demo PRIVATE ../../include/)
