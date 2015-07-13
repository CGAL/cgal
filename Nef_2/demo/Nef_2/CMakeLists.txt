# This is the CMake script for compiling a CGAL application.

project( Nef_2 ) 

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with Qt5: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND )

  include_directories(BEFORE ../../include)
  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/nef_2.cpp" nef_2.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Qt_widget_toolbar.h" Qt_widget_toolbar.moc )


  # The executable itself.
  add_executable  ( nef_2 
                    nef_2.cpp nef_2.moc
		    Qt_widget_toolbar.cpp Qt_widget_toolbar.moc )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS nef_2 )

  target_link_libraries( nef_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
