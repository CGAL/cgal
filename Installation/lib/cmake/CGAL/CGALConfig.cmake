# 
# This file is the CGALConfig.cmake for a pure header-only CGAL installion
#

if(CGALConfig_included)
  return()
endif()

set(CGALConfig_included TRUE)

get_filename_component(CGAL_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

set(CGAL_HEADER_ONLY TRUE)

# Save the current source directory. That variable can be changed by
# a `CMakeLists.txt`, for `CMakeLists.txt` files that are created in
# the binary directory.
set(CGAL_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

function(_detect_branch_build VAR_NAME)
  if(IS_DIRECTORY ${CGAL_CONFIG_DIR}/../../../../Installation)
    set(${VAR_NAME} TRUE PARENT_SCOPE)
  else()
    set(${VAR_NAME} FALSE PARENT_SCOPE)
  endif()
endfunction()

_detect_branch_build(BRANCH_BUILD)
if(BRANCH_BUILD)
  set(CGAL_ROOT ${CGAL_CONFIG_DIR}/../../../..)
  set(CGAL_INSTALLATION_PACKAGE_DIR ${CGAL_ROOT}/Installation)
  set(CGAL_GRAPHICSVIEW_PACKAGE_DIR ${CGAL_ROOT}/GraphicsView)
  file(GLOB include_dirs ${CGAL_ROOT}/*/include)
#  message("inc_dirs: ${include_dirs}")
  foreach(inc_dir ${include_dirs})
    if(IS_DIRECTORY ${inc_dir})
      list(APPEND CGAL_INCLUDE_DIRS ${inc_dir})
      if(EXISTS ${inc_dir}/CGAL/config.h)
	set(CGAL_FOUND TRUE)
      endif()
    endif()
  endforeach()
else()
  set(CGAL_ROOT ${CGAL_CONFIG_DIR}/../../..)
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

set(CGAL_MODULES_DIR ${CGAL_CONFIG_DIR}/../../../cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${CGAL_MODULES_DIR})

include( ${CGAL_MODULES_DIR}/CGAL_SCM.cmake )
CGAL_detect_git(${CGAL_CONFIG_DIR}/../../../..)

include( ${CGAL_MODULES_DIR}/CGAL_SetupCGALDependencies.cmake )

#
# Define the CGAL targets and theirs CGAL:: aliases
#
foreach(cgal_lib CGAL CGAL_Core CGAL_ImageIO CGAL_Qt5)
  set(${cgal_lib}_FOUND TRUE)
  add_library(${cgal_lib} INTERFACE)
  add_library(CGAL::${cgal_lib} ALIAS ${cgal_lib})
  if(NOT ${cgal_lib} STREQUAL CGAL_Core)
    include( ${CGAL_MODULES_DIR}/CGAL_Setup${cgal_lib}Dependencies.cmake )
  endif()
endforeach()
target_compile_definitions(CGAL INTERFACE CGAL_HEADER_ONLY=1)
target_include_directories(CGAL INTERFACE ${CGAL_INCLUDE_DIRS})
CGAL_setup_CGAL_dependencies(CGAL INTERFACE)
CGAL_setup_CGAL_Qt5_dependencies(CGAL_Qt5 INTERFACE)
target_link_libraries( CGAL_Core INTERFACE CGAL::CGAL )
CGAL_setup_CGAL_ImageIO_dependencies(CGAL_ImageIO INTERFACE)

#
#
#

include(Use_CGAL_Qt5_headers)

message("CGAL_FOUND: ${CGAL_FOUND}")
message("CGAL_Qt5_FOUND: ${CGAL_Qt5_FOUND}")
message("Qt5_FOUND: ${Qt5_FOUND}")

include(${CGAL_MODULES_DIR}/CGAL_CreateSingleSourceCGALProgram.cmake)
include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)
