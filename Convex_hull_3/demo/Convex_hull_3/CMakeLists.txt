# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Convex_hull_3_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(${CGAL_USE_FILE})
  include(Qt3Macros)

  QT3_AUTOMOC(  incremental_hull_3_demo.cpp quickhull_3_demo.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )

  add_executable  (quickhull_3_demo  incremental_hull_3_demo.cpp quickhull_3_demo.cpp)
  add_dependencies(quickhull_3_demo CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(quickhull_3_demo CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(quickhull_3_demo ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(quickhull_3_demo ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
