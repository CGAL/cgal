# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Snap_rounding_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(${CGAL_USE_FILE})
  include(Qt3Macros)

  QT3_AUTOMOC(  snap_rounding_2.cpp snap_rounding_2_toolbar.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )

  add_executable  (snap_rounding_2  snap_rounding_2.cpp snap_rounding_2_toolbar.cpp)
  add_dependencies(snap_rounding_2 CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(snap_rounding_2 CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(snap_rounding_2 ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(snap_rounding_2 ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
