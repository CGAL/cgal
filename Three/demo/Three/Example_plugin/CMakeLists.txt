project( Example_plugin )
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Compatibility with CMake 3.0
if(POLICY CMP0042)
  # Do not enable the use of MACOSX_RPATH
  # http://www.cmake.org/cmake/help/v3.0/policy/CMP0042.html
  cmake_policy(SET CMP0042 OLD)
endif()
#Find CGAL
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )
# Find Qt5 itself
find_package(Qt5
             QUIET
             COMPONENTS            OpenGL Script Svg Xml
             OPTIONAL_COMPONENTS   ScriptTools)
             
if(Qt5_FOUND AND CGAL_FOUND)  
  find_package(CGAL_polyhedron_demo)

include( ${CGAL_POLYHEDRON_DEMO_USE_FILE} )

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
           

  polyhedron_demo_plugin(example_plugin Example_plugin)
  
  qt5_wrap_ui( basicUI_FILES Basic_dialog.ui  )
  polyhedron_demo_plugin(basic_plugin Basic_plugin ${basicUI_FILES})

  qt5_wrap_ui( dockUI_FILES Basic_dock_widget.ui )
  polyhedron_demo_plugin(dock_widget_plugin Dock_widget_plugin ${dockUI_FILES})

  polyhedron_demo_plugin(basic_item_plugin Basic_item_plugin)

#  links the library containing the scene_plane_item with the plugin
  target_link_libraries(basic_item_plugin Polyhedron_scene_basic_objects)

endif()
