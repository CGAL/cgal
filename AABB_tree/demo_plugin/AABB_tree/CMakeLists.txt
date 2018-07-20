project( Do_trees_intersect_plugin)
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 3.1)

#Find CGAL
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )
# Find Qt5 itself
find_package(Qt5
             QUIET
             COMPONENTS            OpenGL Script Svg Xml
             OPTIONAL_COMPONENTS   ScriptTools)
find_package(CGAL_polyhedron_demo
HINTS "${CGAL_DIR}" "${CGAL_DIR}/Polyhedron/demo/build-ho"
)
include( ${CGAL_POLYHEDRON_DEMO_USE_FILE} )
#macro add_item
macro(add_item item_name)
  add_library(${item_name} SHARED ${ARGN})
  target_link_libraries(${item_name}
    PUBLIC Polyhedron_demo_framework ${CGAL_LIBRARIES} ${Boost_LIBRARIES}
    Qt5::OpenGL Qt5::Gui Qt5::Script Qt5::Widgets)
  cgal_add_compilation_test(${item_name})
endmacro(add_item)

add_item(scene_movable_sm_item Scene_movable_sm_item.cpp)



polyhedron_demo_plugin(do_trees_intersect_plugin Do_trees_intersect_plugin)
target_link_libraries(do_trees_intersect_plugin PUBLIC Polyhedron_scene_surface_mesh_item scene_movable_sm_item)

