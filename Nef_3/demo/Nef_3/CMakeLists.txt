# This is the CMake script for compiling a CGAL application.

project (Nef_3_Demo)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL QUIET COMPONENTS Qt3)
include(${CGAL_USE_FILE})

find_package(Qt3-patched QUIET )
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )

  include(Qt3Macros-patched)
  include(CreateSingleSourceCGALQt3Program)

  include_directories (BEFORE include)

  create_single_source_cgal_qt3_program( nef_3_extended.cpp    )
  create_single_source_cgal_qt3_program( nef_3_filtered.cpp    )
  create_single_source_cgal_qt3_program( nef_3_homogeneous.cpp )
  create_single_source_cgal_qt3_program( visual_hull.cpp       )
  create_single_source_cgal_qt3_program( visualization_SM.cpp  )
  create_single_source_cgal_qt3_program( visualization_SNC.cpp )

else( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )
