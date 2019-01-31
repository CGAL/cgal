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

find_package(Qt5 QUIET COMPONENTS Widgets)

find_package( Boost REQUIRED COMPONENTS thread )
include_directories(${Boost_INCLUDE_DIRS})
include_directories (BEFORE include)

if(Boost_FOUND AND CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND CGAL_Core_FOUND)
  # ui files, created with Qt Designer
  qt5_wrap_ui( UIS HDT2.ui )
  
  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( RESOURCE_FILES resources/File.qrc resources/Input.qrc resources/Delaunay_triangulation_2.qrc )
  
  # cpp files
  message("RESOURCES :  ${RESOURCE_FILES}")
  add_executable ( HDT2_demo HDT2.cpp ${RESOURCE_FILES} ${UIS})
  target_link_libraries ( HDT2_demo CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core)
  qt5_use_modules( HDT2_demo Widgets )
else()
  message(STATUS "NOTICE: This demo requires CGAL, Boost and Qt5 and will not be compiled.")
endif()
