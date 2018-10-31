cmake_minimum_required(VERSION 2.8)

project(Boolean_set_operations_2_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL COMPONENTS Qt5 Core )

include(${CGAL_USE_FILE})

if (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND CGAL_Core_FOUND)

  set(QT_USE_QTMAIN TRUE)

  find_package(Qt5 QUIET COMPONENTS Widgets Script Svg OpenGL)

  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include )
  include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
  # include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  # include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/include)

  qt5_add_resources(UI_RESOURCES UiFiles/Boolean_operations_2.qrc)

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
    ${UI_HEADERS}
    ${UI_RESOURCES})

  # Link with Qt libraries
  target_link_libraries(Boolean_set_operations_2 ${QT_LIBRARIES})

  # Link with CGAL
  #target_link_libraries(Boolean_set_operations_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

  qt5_use_modules(Boolean_set_operations_2 Widgets Script Svg OpenGL)

  #set_target_properties(Boolean_set_operations_2 PROPERTIES LINK_FLAGS "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")

else()

  message(STATUS "NOTICE: This demo requires CGAL, CGAL-Core and Qt5, and will not be compiled.")

endif()
