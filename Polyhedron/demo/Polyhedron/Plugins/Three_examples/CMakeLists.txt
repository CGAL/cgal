cmake_minimum_required(VERSION 3.1...3.15)
project(Three_examples)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
#Find CGAL
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5 ImageIO)
# Find Qt5 itself
find_package(
  Qt5 QUIET
  COMPONENTS OpenGL Script Svg Xml
  OPTIONAL_COMPONENTS ScriptTools WebSockets)

if(RUNNING_CGAL_AUTO_TEST)
  if(Qt5_FOUND)
    include(${CGAL_USE_FILE})
  endif()
  polyhedron_demo_plugin(example_plugin Example_plugin)

  qt5_wrap_ui(basicUI_FILES Basic_dialog.ui)
  polyhedron_demo_plugin(basic_plugin Basic_plugin ${basicUI_FILES})

  qt5_wrap_ui(dockUI_FILES Basic_dock_widget.ui)
  polyhedron_demo_plugin(dock_widget_plugin Dock_widget_plugin ${dockUI_FILES})

  polyhedron_demo_plugin(basic_item_plugin Basic_item_plugin)

  #  links the library containing the scene_plane_item with the plugin
  target_link_libraries(basic_item_plugin PUBLIC scene_basic_objects)
endif()
