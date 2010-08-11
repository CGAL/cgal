# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Boolean_operations_2_demo ) 

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt4 Core )

include(${CGAL_USE_FILE})

set( QT_USE_QTMAIN   TRUE )

find_package(Qt4)

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include )

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND CGAL_Core_FOUND )

  include(${QT_USE_FILE})

  # UI files (Qt Designer files)
  qt4_wrap_ui( CDT_UI_FILES boolean_operations_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( CDT_RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt4_generate_moc( "boolean_operations_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/boolean_operations_2.moc" )

  # The executable itself.
  add_executable  ( boolean_operations_2 ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.cpp boolean_operations_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

  # Link with Qt libraries
  target_link_libraries( boolean_operations_2 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( boolean_operations_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
