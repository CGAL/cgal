# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Periodic_4_hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories( ${CMAKE_BINARY_DIR} )

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON) 

cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL QUIET COMPONENTS Core Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Widgets)
if(CGAL_FOUND AND CGAL_Core_FOUND AND Qt5_FOUND AND CGAL_Qt5_FOUND)

  find_package( Boost REQUIRED )
  include_directories(${Boost_INCLUDE_DIRS})

  include_directories (BEFORE include )

  # ui files, created with Qt Designer
  qt5_wrap_ui( UIS P4HDT2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( RESOURCE_FILES Main_resources.qrc )

  # cpp files

  add_executable ( P4HDT2_demo
  P4HDT2.cpp ${RESOURCE_FILES} ${UIS})

  #add_executable ( Periodic_4_hyperbolic_billiards_demo
  #Periodic_4_hyperbolic_billiards_demo.cpp ${RESOURCE_FILES} )

  qt5_use_modules( P4HDT2_demo Widgets )
  #qt5_use_modules( Periodic_4_hyperbolic_billiards_demo Widgets )

  target_link_libraries( P4HDT2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

  #target_link_libraries( Periodic_4_hyperbolic_billiards_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
else()
  message(STATUS "NOTICE: This demo requires Qt5 and will not be compiled.")
endif()
