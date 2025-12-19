# This is the CMake script for compiling the CGAL Mesh_3 demo implicit functions.

cmake_minimum_required(VERSION 3.12...3.31)

include(CGALlab_macros)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include directory of demo includes
include_directories(BEFORE ${Mesh_3_implicit_functions_BINARY_DIR} ../include)

# Find CGAL and CGAL Qt6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Core)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  ###########
  # PLUGINS #
  ###########
  include(AddFileDependencies)
  remove_definitions(-DQT_STATICPLUGIN)

  cgal_lab_plugin(p_sphere_function_plugin Sphere_implicit_function
                         KEYWORDS Mesh_3)
  cgal_lab_plugin(p_tanglecube_function_plugin
                         Tanglecube_implicit_function KEYWORDS Mesh_3)
  cgal_lab_plugin(p_klein_function_plugin Klein_implicit_function
                         KEYWORDS Mesh_3)

else(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt6 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(MESH_3_MISSING_DEPS "Qt6, ${MESH_3_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${MESH_3_MISSING_DEPS} and will not be compiled.")

endif(CGAL_Qt6_FOUND AND Qt6_FOUND)
