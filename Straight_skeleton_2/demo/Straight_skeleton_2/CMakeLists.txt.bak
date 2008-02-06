project (Straight_skeleton_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_BINARY_DIR}/cmake/modules" )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules"  )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)
find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(UseCGAL)
  include(Qt3Macros)
  include(CreateSingleSourceCGALQt3Program)
  
  
  # Make sure the compiler can find local include files
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    include_directories (BEFORE include)
  endif()
  
  create_single_source_cgal_qt3_program( straight_skeleton_2.cpp straight_skeleton_2_toolbar.cpp straight_skeleton_2_toolbar_layers.cpp)
  
endif()
