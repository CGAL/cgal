# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Arrangement_on_surface_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(${CGAL_USE_FILE})
  include(Qt3Macros)

  QT3_AUTOMOC(  arrangement_2.cpp demo_tab.cpp forms.cpp MyWindow_files.cpp MyWindow_operations.cpp MyWindow_overlay.cpp qt_layer.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )

  add_executable  (arrangement_2  arrangement_2.cpp demo_tab.cpp forms.cpp MyWindow_files.cpp MyWindow_operations.cpp MyWindow_overlay.cpp qt_layer.cpp)
  add_dependencies(arrangement_2 CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(arrangement_2 CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(arrangement_2 ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(arrangement_2 ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
