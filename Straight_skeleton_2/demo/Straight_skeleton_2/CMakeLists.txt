# This is the CMake script for compiling a CGAL application.

project( straight_skeleton_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt3 Core )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND )

  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/straight_skeleton_2_toolbar.h" straight_skeleton_2_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/straight_skeleton_2_toolbar_layers.h" straight_skeleton_2_toolbar_layers.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/straight_skeleton_2.cpp" straight_skeleton_2.moc )


  # The executable itself.
  add_executable  ( straight_skeleton_2 
    straight_skeleton_2.cpp straight_skeleton_2.moc
    straight_skeleton_2_toolbar.cpp straight_skeleton_2_toolbar.moc 
    straight_skeleton_2_toolbar_layers.cpp straight_skeleton_2_toolbar_layers.moc
    )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS straight_skeleton_2 )

  target_link_libraries( straight_skeleton_2  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT3_LIBRARIES} )
  
message( STATUS "QT3_LIBRARIES=${QT3_LIBRARIES}" )   
message( STATUS "CGAL_LIBRARIES=${CGAL_LIBRARIES}" )   
message( STATUS "CGAL_3RD_PARTY_LIBRARIES=${CGAL_3RD_PARTY_LIBRARIES}" )   

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
