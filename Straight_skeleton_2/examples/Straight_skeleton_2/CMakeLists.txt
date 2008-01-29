project (Straight_skeleton_2_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_BINARY_DIR}/cmake/modules" )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)
  
  # Make sure the compiler can find local include files
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    include_directories (BEFORE include)
  endif()
  
  create_single_source_cgal_program( Create_offset_polygons_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_from_polygon_with_holes_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_with_holes_2.cpp )
  create_single_source_cgal_program( Create_straight_skeleton_2.cpp )
  create_single_source_cgal_program( Create_straight_skeleton_from_polygon_with_holes_2.cpp )
  create_single_source_cgal_program( Low_level_API.cpp )
  
endif()



