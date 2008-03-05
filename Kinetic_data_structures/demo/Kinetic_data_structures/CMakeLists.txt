# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kinetic_data_structures_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(${CGAL_USE_FILE})
  include(Qt3Macros)

  QT3_AUTOMOC(  Delaunay_triangulation_2.cpp Delaunay_triangulation_3.cpp Delaunay_triangulation_stable_subset_2.cpp generate_data.cpp gui_2.cpp gui_3.cpp regular_triangulation_3.cpp SoQt_examiner_viewer.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )
  include_directories (BEFORE include)

  add_executable  (regular_triangulation_3  Delaunay_triangulation_2.cpp Delaunay_triangulation_3.cpp Delaunay_triangulation_stable_subset_2.cpp generate_data.cpp gui_2.cpp gui_3.cpp regular_triangulation_3.cpp SoQt_examiner_viewer.cpp)
  add_dependencies(regular_triangulation_3 CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(regular_triangulation_3 CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(regular_triangulation_3 ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(regular_triangulation_3 ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
