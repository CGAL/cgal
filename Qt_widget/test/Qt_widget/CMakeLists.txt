# This is the CMake script for compiling a CGAL application.

project( Qt_widget_test ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)

find_package(CGAL QUIET COMPONENTS Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched QUIET)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND )


  # The executable itself.
  add_executable  ( test_qt_version
                    test_qt_version.cpp )

  target_link_libraries( test_qt_version ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
