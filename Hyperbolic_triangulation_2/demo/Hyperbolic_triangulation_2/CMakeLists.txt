# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON) 

cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL QUIET COMPONENTS Core Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Widgets)

find_package( Boost REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

include_directories (BEFORE include)

# ui files, created with Qt Designer
qt5_wrap_ui( UIS Hyperbolic_triangulation_2_demo.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( RESOURCE_FILES resources/File.qrc resources/Input.qrc resources/Delaunay_triangulation_2.qrc )

# cpp files

add_executable ( Hyperbolic_Delaunay_triangulation_2_demo Hyperbolic_Delaunay_triangulation_2_demo.cpp ${RESOURCE_FILES} ${UIS})

qt5_use_modules( Hyperbolic_Delaunay_triangulation_2_demo Widgets )

