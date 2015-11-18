project( Example_plugin )
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
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
  find_package(CGAL_polyhedron_demo
  HINTS "${CGAL_DIR}" "${CGAL_DIR}/Polyhedron/demo/Polyhedron-build"
  )

include( ${CGAL_POLYHEDRON_DEMO_USE_FILE} )
# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

           
  #include( ${CMAKE_CURRENT_SOURCE_DIR}/../../Polyhedron/demo/Polyhedron/polyhedron_demo_macros.cmake )
  qt5_wrap_ui( exampleUI_FILES Polyhedron_demo_example_plugin.ui dock_example.ui)
  qt5_wrap_ui( example_dockUI_FILES dock_example.ui)
  polyhedron_demo_plugin(example_plugin Polyhedron_demo_example_plugin ${exampleUI_FILES} ${example_dockUI_FILES})
  target_link_libraries(example_plugin scene_polyhedron_item)  
endif()