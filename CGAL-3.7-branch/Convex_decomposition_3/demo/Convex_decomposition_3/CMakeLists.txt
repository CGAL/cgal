# This is the CMake script for compiling a CGAL application.

project( Convex_decomposition_3 ) 

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


  # The executable itself.
  add_executable  ( convex_decomposition_3
                    convex_decomposition_3.cpp )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS convex_decomposition_3 )

  target_link_libraries( convex_decomposition_3 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
