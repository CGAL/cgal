# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Interpolation_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(UseCGAL)
  include(Qt3Macros)

  QT3_AUTOMOC(  interpolation_2.cpp surface_voronoi.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )

  add_executable  (surface_voronoi  interpolation_2.cpp surface_voronoi.cpp)
  add_dependencies(surface_voronoi CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(surface_voronoi CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(surface_voronoi ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(surface_voronoi ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
