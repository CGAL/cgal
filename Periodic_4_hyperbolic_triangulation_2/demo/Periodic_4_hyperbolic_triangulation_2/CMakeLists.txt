# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (My_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON) 
set(CMAKE_AUTOUIC ON)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL QUIET COMPONENTS Core Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Widgets)

find_package( Boost REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

include_directories (BEFORE ../../include include )

# ui files, created with Qt Designer
qt5_wrap_ui( uis Periodic_4_hyperbolic_Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( RESOURCE_FILES Main_resources.qrc )


# cpp files

#add_executable ( Playground playground.cpp )
#add_executable ( Triangulation triangulation.cpp )
add_executable ( Demo Periodic_4_hyperbolic_Delaunay_triangulation_2_demo.cpp ${RESOURCE_FILES} )
add_executable ( Overlay_demo Periodic_4_hyperbolic_Delaunay_triangulation_2_overlay_demo.cpp ${RESOURCE_FILES} )
qt5_use_modules( Demo Widgets )
qt5_use_modules( Overlay_demo Widgets )
#target_link_libraries( Playground ${CGAL_LIBRARIES} )
#target_link_libraries( Triangulation ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
target_link_libraries( Demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
target_link_libraries( Overlay_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

