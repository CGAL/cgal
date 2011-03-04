# This is the CMake script for compiling a CGAL application.

project( Stream_lines_2 ) 

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

find_package(OpenGL)

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND AND OPENGL_FOUND)

  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/streamlines.cpp" streamlines.moc )

  include_directories(${OPENGL_INCLUDE_DIR})

  # The executable itself.
  add_executable  ( streamlines
                    streamlines.cpp streamlines.moc )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS streamlines )

  target_link_libraries( streamlines ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES} ${OPENGL_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL, Qt3, and OpenGL, and will not be compiled.")

endif()
