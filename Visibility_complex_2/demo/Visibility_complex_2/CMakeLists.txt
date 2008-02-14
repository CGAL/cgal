# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Visibility_complex_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

find_package(Qt3)

if ( CGAL_FOUND AND QT_FOUND )

  include(UseCGAL)
  include(Qt3Macros)

  QT3_AUTOMOC(  arrangements.cpp circle_compute_bit.cpp circle_compute_vc.cpp circle_display.cpp circle_input.cpp ellipse_compute_bit.cpp ellipse_compute_vc.cpp ellipse_display.cpp ellipse_input.cpp point_compute_bit.cpp point_compute_sp.cpp point_compute_vc.cpp point_display.cpp point_input.cpp polygon_compute_bit.cpp polygon_compute_sp.cpp polygon_compute_vc.cpp polygon_display.cpp polygon_input.cpp Qt_widget_styled_layer.cpp Qt_widget_style_editor.cpp segment_compute_bit.cpp segment_compute_sp.cpp segment_compute_vc.cpp segment_display.cpp segment_input.cpp Show_lines_base.cpp Show_points.cpp Show_pseudo_triangulation.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )
  include_directories (BEFORE include)

  add_executable  (segment_input  arrangements.cpp circle_compute_bit.cpp circle_compute_vc.cpp circle_display.cpp circle_input.cpp ellipse_compute_bit.cpp ellipse_compute_vc.cpp ellipse_display.cpp ellipse_input.cpp point_compute_bit.cpp point_compute_sp.cpp point_compute_vc.cpp point_display.cpp point_input.cpp polygon_compute_bit.cpp polygon_compute_sp.cpp polygon_compute_vc.cpp polygon_display.cpp polygon_input.cpp Qt_widget_styled_layer.cpp Qt_widget_style_editor.cpp segment_compute_bit.cpp segment_compute_sp.cpp segment_compute_vc.cpp segment_display.cpp segment_input.cpp Show_lines_base.cpp Show_points.cpp Show_pseudo_triangulation.cpp)
  add_dependencies(segment_input CGAL CGAL_QT)

  if ( CGAL_USE_CGAL_CORE )
    add_dependencies(segment_input CGAL_CORE)
  endif()

  # Link the executable to CGAL and third-party libraries
  if ( AUTO_LINK_ENABLED )
    target_link_libraries(segment_input ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
  else()
    target_link_libraries(segment_input ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
  endif()

endif()
