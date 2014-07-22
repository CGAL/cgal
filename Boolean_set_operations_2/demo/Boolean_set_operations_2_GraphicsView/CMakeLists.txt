# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Boolean_operations_2_demo ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set( QT_USE_QTMAIN   TRUE )

#New for Qt5 version !
option (USE_QT5 "Use Qt5 version instead of Qt4" ON)

if(USE_QT5)

	find_package(CGAL COMPONENTS Qt5 Core )

	include(${CGAL_USE_FILE})
	
	find_package(Qt5)

else (USE_QT5)

	
	find_package(CGAL COMPONENTS Qt4 Core )
	include(${CGAL_USE_FILE})

	find_package(Qt4)

	if(CGAL_Qt4_FOUND AND QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
	endif()

endif (USE_QT5)

include(QtChoice)

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include )

if ( CGAL_FOUND AND (QT4 OR QT5) AND CGAL_Core_FOUND )

  #include(${QT_USE_FILE})

  # UI files (Qt Designer files)
  qt_wrap_ui( CDT_UI_FILES boolean_operations_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt_add_resources ( CDT_RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt_generate_moc( "boolean_operations_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/boolean_operations_2.moc" )

  # The executable itself.
  add_executable  ( boolean_operations_2 ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.cpp boolean_operations_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

  # Link with Qt libraries
  target_link_libraries( boolean_operations_2 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( boolean_operations_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else( CGAL_FOUND AND (QT4 OR QT5) AND CGAL_Core_FOUND )

  message(STATUS "NOTICE: This demo requires CGAL and Qt ( version 4 or 5), and will not be compiled.")

endif( CGAL_FOUND AND (QT4 OR QT5) AND CGAL_Core_FOUND )
