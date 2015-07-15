# This is the CMake script for compiling the CGAL Mesh_3 demo implicit functions.

project( Mesh_3_implicit_functions )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include directory of demo includes
include_directories( BEFORE ${Mesh_3_implicit_functions_BINARY_DIR} ../include )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL ScriptTools )


if(CGAL_Qt5_FOUND AND Qt5_FOUND)
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
    qt5_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    qt5_use_modules(${plugin_name} Widgets Script OpenGL Gui Xml)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(p_sphere_function_plugin Sphere_implicit_function)
  polyhedron_demo_plugin(p_tanglecube_function_plugin Tanglecube_implicit_function)
  polyhedron_demo_plugin(p_klein_function_plugin Klein_implicit_function)

else (CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt5 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(MESH_3_MISSING_DEPS "Qt5, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt5_FOUND AND Qt5_FOUND)
