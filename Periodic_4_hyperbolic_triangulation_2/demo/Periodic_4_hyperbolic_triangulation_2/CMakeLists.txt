# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Periodic_4_hyperbolic_triangulation_2_demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTOUIC ON)  

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Widgets)

include_directories (BEFORE ../../include include ../../../Hyperbolic_triangulation_2/include)

# ui files, created with Qt Designer
qt5_wrap_ui( uis Periodic_4_hyperbolic_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( RESOURCE_FILES resources/Periodic_4_hyperbolic_triangulation_2.qrc )


# cpp files
add_executable ( Periodic_4_hyperbolic_triangulation_2_demo Periodic_4_hyperbolic_triangulation_2_demo.cpp )
qt5_use_modules( Periodic_4_hyperbolic_triangulation_2_demo Widgets)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_4_hyperbolic_triangulation_2_demo )
target_link_libraries( Periodic_4_hyperbolic_triangulation_2_demo ${CGAL_LIBRARIES} ${QT_LIBRARIES} ${CGAL_QT_LIBRARIES} )


