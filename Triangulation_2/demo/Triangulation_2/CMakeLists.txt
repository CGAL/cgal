# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(UseCGAL)
  include(Qt3Macros)

  QT3_AUTOMOC(  constrained.cpp constrained_delaunay_triangulation_2.cpp constrained_delaunay_triangulation_2_toolbar.cpp constrained_delaunay_triangulation_2_toolbar_layers.cpp delaunay_triangulation_2.cpp delaunay_triangulation_2_toolbar.cpp delaunay_triangulation_2_toolbar_layers.cpp regular_triangulation_2.cpp regular_triangulation_2_toolbar.cpp regular_triangulation_2_toolbar_layers.cpp triangulation_2_edit_vertex.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )

  add_executable  (regular_triangulation_2  constrained.cpp constrained_delaunay_triangulation_2.cpp constrained_delaunay_triangulation_2_toolbar.cpp constrained_delaunay_triangulation_2_toolbar_layers.cpp delaunay_triangulation_2.cpp delaunay_triangulation_2_toolbar.cpp delaunay_triangulation_2_toolbar_layers.cpp regular_triangulation_2.cpp regular_triangulation_2_toolbar.cpp regular_triangulation_2_toolbar_layers.cpp triangulation_2_edit_vertex.cpp)
  add_dependencies(regular_triangulation_2 CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(regular_triangulation_2 CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(regular_triangulation_2 ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(regular_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
