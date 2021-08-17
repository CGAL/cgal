cmake_minimum_required(VERSION 2.8)

project(Boolean_set_operations_2_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()
if(POLICY CMP0072)
  # About the use of OpenGL
  cmake_policy(SET CMP0072 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
if (has_cpp11 LESS 0)
  message(STATUS "NOTICE: This demo requires a C++11 compiler and will not be compiled.")
  return()
endif()

# Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL COMPONENTS Core Qt5 )

include(${CGAL_USE_FILE})

if (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

  set(QT_USE_QTMAIN TRUE)

  find_package(Qt5 QUIET COMPONENTS Widgets Script Svg OpenGL)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include )

  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES Boolean_operations_2.qrc)

  qt5_wrap_ui(UI_HEADERS
    UiFiles/Boolean_set_operations_2.ui
    UiFiles/ChoiseDialogPolygon.ui)

  add_executable(Boolean_set_operations_2
    Boolean_set_operations_2.cpp
    # MainWindow.cpp
    # include/MainWindow.h
    include/Typedefs.h
    include/QT5/Boundary_pieces_graphics_item.h
    include/QT5/Circular_polygons.h
    include/QT5/Graphics_view_circular_polygon_input.h
    include/QT5/Piecewise_boundary_graphics_item.h
    include/QT5/Piecewise_graphics_item_base.h
    include/QT5/Piecewise_region_graphics_item.h
    include/QT5/Piecewise_set_graphics_item.h
    include/QT5/Graphics_view_linear_polygon_input.h
    include/QT5/Linear_polygons.h
    include/QT5/General_polygon_2.h
    include/QT5/General_polygon_set_2.h
    include/QT5/General_polygon_set_on_surface_2.h
    include/QT5/Gps_circle_segment_traits_2.h
    include/QT5/Gps_segment_traits_2.h
    include/QT5/Gps_traits_2.h
    include/QT5/connect_holes.h
    include/QT5/Polygon_set_2.h
    ${UI_HEADERS}
    ${CGAL_Qt5_RESOURCE_FILES})

  # Link with CGAL and Qt libraries
  target_link_libraries(Boolean_set_operations_2
    PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets Qt5::Script Qt5::Svg Qt5::OpenGL)

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5 (Widges, Script, Svg, OpenGL), and will not be compiled.")

endif()
