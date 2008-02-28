# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project (Nef_3_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(${CGAL_USE_FILE})
  include(Qt3Macros)
  include(CreateSingleSourceCGALQt3Program)

  include_directories (BEFORE include)

  create_single_source_cgal_qt3_program( nef_3_extended.cpp    )
  create_single_source_cgal_qt3_program( nef_3_filtered.cpp    )
  create_single_source_cgal_qt3_program( nef_3_homogeneous.cpp )
  create_single_source_cgal_qt3_program( visual_hull.cpp       )
  create_single_source_cgal_qt3_program( visualization_SM.cpp  )
  create_single_source_cgal_qt3_program( visualization_SNC.cpp )

else()

  message( STATUS "This demo requires CGAL and Qt3" )
  
endif()
