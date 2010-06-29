# This is the CMake script for compiling the CGAL Mesh_3 demo implicit functions.

project( Mesh_3_implicit_functions )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Include directory of demo includes
include_directories( BEFORE ${Mesh_3_implicit_functions_BINARY_DIR} ../include )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)


if(CGAL_Qt4_FOUND AND QT4_FOUND)
  include(${QT_USE_FILE})

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  ###########
  # PLUGINS #
  ###########
  include(AddFileDependencies)
  remove_definitions(-DQT_STATICPLUGIN)

  # polyhedron_demo_plugin is a macro copied from demo/Polyhedron/CMakeLists.txt
  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      set(other_sources ${ARGN})
      set(option "")
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(sphere_function_plugin Sphere_implicit_function)
  polyhedron_demo_plugin(tanglecube_function_plugin Tanglecube_implicit_function)
  polyhedron_demo_plugin(klein_function_plugin Klein_implicit_function)

else (CGAL_Qt4_FOUND AND QT4_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt4 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(MESH_3_MISSING_DEPS "Qt4, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND)
