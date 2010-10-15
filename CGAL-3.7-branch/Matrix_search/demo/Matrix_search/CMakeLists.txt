# This is the CMake script for compiling a CGAL application.

project( Matrix_search_Demo ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND)

  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Qt_widget_toolbar.h" Qt_widget_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Qt_widget_move_list_point.h" Qt_widget_move_list_point.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/max_k-gon.cpp" max_k-gon.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/rectangular_p_center_2_demo.cpp" rectangular_p_center_2_demo.moc )


  # The executable itself.
  add_executable  ( max_k-gon
    max_k-gon.cpp max_k-gon.moc
    Qt_widget_toolbar.cpp Qt_widget_toolbar.moc
    Qt_widget_move_list_point.cpp Qt_widget_move_list_point.moc
    )

  # The executable itself.
  add_executable  ( rectangular_p_center_2_demo
    rectangular_p_center_2_demo.cpp rectangular_p_center_2_demo.moc
    Qt_widget_toolbar.cpp Qt_widget_toolbar.moc
    Qt_widget_move_list_point.cpp Qt_widget_move_list_point.moc
    )

  target_link_libraries( max_k-gon  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  target_link_libraries( rectangular_p_center_2_demo  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS max_k-gon )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS rectangular_p_center_2_demo )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
