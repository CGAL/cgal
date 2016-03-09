# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Boolean_operations_2_demo ) 

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5 Core )

include(${CGAL_USE_FILE})

set( QT_USE_QTMAIN   TRUE )

find_package(Qt5 QUIET COMPONENTS Widgets Script Svg)

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include )

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND CGAL_Core_FOUND )

  # UI files (Qt Designer files)
  qt5_wrap_ui( CDT_UI_FILES boolean_operations_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt5_generate_moc( "boolean_operations_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/boolean_operations_2.moc" )

  # The executable itself.
  add_executable  ( boolean_operations_2 ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.cpp boolean_operations_2.moc ${CGAL_Qt5_MOC_FILES} ${CDT_UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} )
  qt5_use_modules(boolean_operations_2 Widgets Script Svg)
  # Link with Qt libraries
  target_link_libraries( boolean_operations_2 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( boolean_operations_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
