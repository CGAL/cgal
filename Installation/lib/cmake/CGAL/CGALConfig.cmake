# 
# This file is the CGALConfig.cmake for a header-only CGAL installation
#

# For UseCGAL.cmake
set( CGAL_REQUESTED_COMPONENTS ${CGAL_FIND_COMPONENTS} )

set(CGAL_LIBRARIES CGAL)

if(CGAL_BUILDING_LIBS)
  foreach(comp ${CGAL_FIND_COMPONENTS})
    if(CGAL_${comp}_FOUND)
      list(APPEND CGAL_LIBRARIES CGAL_${comp})
    endif()
  endforeach()
  return()
endif()

foreach(comp ${CGAL_FIND_COMPONENTS})
  if(NOT comp MATCHES "Core|ImageIO|Qt5")
    message(FATAL_ERROR "The requested CGAL component ${comp} does not exist!")
  endif()
  list(APPEND CGAL_LIBRARIES CGAL_${comp})
endforeach()

set(CGALConfig_all_targets_are_defined TRUE)
foreach(cgal_lib ${CGAL_LIBRARIES})
  if(NOT WITH_${cgal_lib})
    set(CGALConfig_all_targets_are_defined FALSE)
  endif()
endforeach()
if(CGALConfig_all_targets_are_defined)
  return()
endif()

message(STATUS "Using header-only CGAL")

get_filename_component(CGAL_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

set(CGAL_HEADER_ONLY TRUE)

# Save the current source directory. That variable can be changed by
# a `CMakeLists.txt`, for `CMakeLists.txt` files that are created in
# the binary directory.
set(CGAL_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

function(cgal_detect_branch_build VAR_NAME)
  if(IS_DIRECTORY ${CGAL_CONFIG_DIR}/../../../../Installation/package_info/Installation/)
    set(${VAR_NAME} TRUE PARENT_SCOPE)
  else()
    set(${VAR_NAME} FALSE PARENT_SCOPE)
  endif()
endfunction()

cgal_detect_branch_build(BRANCH_BUILD)
if(BRANCH_BUILD)
  set(CGAL_ROOT ${CGAL_CONFIG_DIR}/../../../..)
  set(CGAL_INSTALLATION_PACKAGE_DIR ${CGAL_ROOT}/Installation)
  set(CGAL_GRAPHICSVIEW_PACKAGE_DIR ${CGAL_ROOT}/GraphicsView)
  set(CGAL_MODULES_DIR ${CGAL_ROOT}/Installation/cmake/modules)
  file(GLOB packages_dirs ${CGAL_ROOT}/*)
#  message("packages_dirs: ${packages_dirs}")
  foreach(package_dir ${packages_dirs})
    set(inc_dir ${package_dir}/include)
    if(IS_DIRECTORY ${inc_dir}
	AND IS_DIRECTORY ${package_dir}/package_info)
      list(APPEND CGAL_INCLUDE_DIRS ${inc_dir})
      if(EXISTS ${inc_dir}/CGAL/config.h)
	set(CGAL_FOUND TRUE)
      endif()
    endif()
  endforeach()
else()
  set(CGAL_ROOT ${CGAL_CONFIG_DIR}/../../..)

  # not BRANCH_BUILD: it can be an installed CGAL, or the tarball layout
  if(EXISTS ${CGAL_CONFIG_DIR}/CGAL_add_test.cmake)
    # installed CGAL
    set(CGAL_MODULES_DIR ${CGAL_CONFIG_DIR})
  else()
    # tarball
    set(CGAL_MODULES_DIR ${CGAL_ROOT}/cmake/modules)
  endif()

  set(CGAL_INSTALLATION_PACKAGE_DIR ${CGAL_ROOT})
  set(CGAL_GRAPHICSVIEW_PACKAGE_DIR ${CGAL_ROOT})
  set(CGAL_INCLUDE_DIRS ${CGAL_ROOT}/include)
  if(EXISTS ${CGAL_ROOT}/include/CGAL/config.h)
    set(CGAL_FOUND TRUE)
  endif()
endif()

if(NOT CGAL_FOUND)
  return()
endif()

list(APPEND CMAKE_MODULE_PATH ${CGAL_MODULES_DIR})

include( ${CGAL_MODULES_DIR}/CGAL_SCM.cmake )
CGAL_detect_git(${CGAL_CONFIG_DIR}/../../../..)

#
# Search for all dependencies
#
foreach(cgal_lib ${CGAL_LIBRARIES})
  include(CGAL_Setup${cgal_lib}Dependencies)
endforeach()

# this include has to be after the loop that includes the
# `CGAL_Setup${cgal_lib}Dependencies` files
include(CGAL_setup_target_dependencies)

#
# Define the CGAL targets and their CGAL:: aliases
#
foreach(cgal_lib ${CGAL_LIBRARIES})
  set(WITH_${cgal_lib} TRUE)
  if(${cgal_lib}_FOUND AND NOT TARGET ${cgal_lib})
    add_library(${cgal_lib} INTERFACE)
    if(NOT TARGET CGAL::${cgal_lib})
      add_library(CGAL::${cgal_lib} ALIAS ${cgal_lib})
    endif()
    if(${cgal_lib} STREQUAL CGAL)
      target_compile_definitions(CGAL INTERFACE CGAL_HEADER_ONLY=1)
    endif()
    CGAL_setup_target_dependencies(${cgal_lib} INTERFACE)
  endif()
endforeach()


#
#
#

include(${CGAL_MODULES_DIR}/CGAL_CreateSingleSourceCGALProgram.cmake)
include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)

# Temporary? Change the CMAKE module path
cgal_setup_module_path()

set(CGAL_USE_FILE ${CGAL_MODULES_DIR}/UseCGAL.cmake)

include("${CGAL_MODULES_DIR}/CGAL_parse_version_h.cmake")
cgal_parse_version_h( "${CGAL_INSTALLATION_PACKAGE_DIR}/include/CGAL/version.h"
  "CGAL_VERSION_NAME"
  "CGAL_MAJOR_VERSION"
  "CGAL_MINOR_VERSION"
  "CGAL_BUGFIX_VERSION"
  "CGAL_BUILD_VERSION")
set(CGAL_VERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}.${CGAL_BUILD_VERSION}")

include("${CGAL_MODULES_DIR}/CGAL_enable_end_of_configuration_hook.cmake")
