#
# This file is the CGALConfig.cmake for a header-only CGAL installation
#

# For UseCGAL.cmake
set( CGAL_REQUESTED_COMPONENTS ${CGAL_FIND_COMPONENTS} )

set(CGAL_LIBRARIES CGAL)

get_filename_component(CGAL_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

function(cgal_detect_branch_build VAR_NAME)
  if(IS_DIRECTORY ${CGAL_CONFIG_DIR}/../../../../Installation/package_info/Installation/)
    set(${VAR_NAME} TRUE PARENT_SCOPE)
  else()
    set(${VAR_NAME} FALSE PARENT_SCOPE)
  endif()
endfunction()

set(CGAL_FOUND FALSE)
cgal_detect_branch_build(BRANCH_BUILD)
if(BRANCH_BUILD)
  set(CGAL_ROOT ${CGAL_CONFIG_DIR})
  get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
  get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
  get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
  get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
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
  include(${CGAL_CONFIG_DIR}/CGALConfig-installation-dirs.cmake OPTIONAL RESULT_VARIABLE _has_installation_dirs)
  if(NOT _has_installation_dirs)
    set(CGAL_ROOT ${CGAL_CONFIG_DIR})
    get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
    if(NOT EXISTS ${CGAL_ROOT}/include/CGAL/config.h)
      get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
    endif()
    if(NOT EXISTS ${CGAL_ROOT}/include/CGAL/config.h)
      get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
    endif()
    if(NOT EXISTS ${CGAL_ROOT}/include/CGAL/config.h)
      get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
    endif()
  endif()
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

#set CGAL_DATA_DIR
if (NOT CGAL_DATA_DIR)
  if(DEFINED ENV{CGAL_DATA_DIR})
    set(CGAL_DATA_DIR $ENV{CGAL_DATA_DIR})
  else()
    if (EXISTS "${CGAL_ROOT}/Data/data")
      set(CGAL_DATA_DIR "${CGAL_ROOT}/Data/data")
    else()
      if (EXISTS "${CGAL_ROOT}/data")
        set(CGAL_DATA_DIR "${CGAL_ROOT}/data")
      else()
        if (EXISTS "${CMAKE_SOURCE_DIR}/../data")
          set(CGAL_DATA_DIR "${CMAKE_SOURCE_DIR}/../data")
        else()
          if (EXISTS "${CMAKE_SOURCE_DIR}/../../data")
              set(CGAL_DATA_DIR "${CMAKE_SOURCE_DIR}/../../data")
          else()
            if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
              message(WARNING "CGAL_DATA_DIR cannot be deduced, set the variable CGAL_DATA_DIR to set the default value of CGAL::data_file_path()")
            endif()
          endif()
        endif()
      endif()
    endif()
  endif()
endif()

if(NOT TARGET CGAL::Data)
  add_library(CGAL::Data INTERFACE IMPORTED GLOBAL)
  if ( NOT "${CGAL_DATA_DIR}" STREQUAL "" )
    set_target_properties(CGAL::Data PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "CGAL_DATA_DIR=\"${CGAL_DATA_DIR}\"")
  endif()
endif()

include(${CGAL_MODULES_DIR}/CGAL_CreateSingleSourceCGALProgram.cmake)
include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)
include(${CGAL_MODULES_DIR}/CGAL_Common.cmake)
include(${CGAL_MODULES_DIR}/CGAL_TweakFindBoost.cmake)
include(${CGAL_MODULES_DIR}/CGAL_enable_end_of_configuration_hook.cmake)

set(CGAL_USE_FILE ${CGAL_MODULES_DIR}/UseCGAL.cmake)

include(${CGAL_CONFIG_DIR}/CGALConfigVersion.cmake)

# Temporary? Change the CMAKE module path
cgal_setup_module_path()

include(${CGAL_MODULES_DIR}/CGAL_target_use_TBB.cmake)

if( CGAL_DEV_MODE OR RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE )
  # Do not use -isystem for CGAL include paths
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
endif()

foreach(comp ${CGAL_FIND_COMPONENTS})
  if(NOT comp MATCHES "Core|ImageIO|Qt5")
    message(FATAL_ERROR "The requested CGAL component ${comp} does not exist!")
  endif()
  if(comp MATCHES "Core" AND CGAL_DISABLE_GMP)
    message("CGAL_Core needs GMP and won't be used.")
  else()
    list(APPEND CGAL_LIBRARIES CGAL_${comp})
  endif()
endforeach()

set(CGALConfig_all_targets_are_defined TRUE)
foreach(cgal_lib ${CGAL_LIBRARIES})
  if(TARGET CGAL::${cgal_lib})
    set(${cgal_lib}_FOUND TRUE)
  else()
    set(CGALConfig_all_targets_are_defined FALSE)
  endif()
endforeach()
if(CGALConfig_all_targets_are_defined)
  return()
endif()

message(STATUS "Using header-only CGAL")

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
    add_library(${cgal_lib} INTERFACE IMPORTED GLOBAL)
    if(NOT TARGET CGAL::${cgal_lib})
      add_library(CGAL::${cgal_lib} ALIAS ${cgal_lib})
    endif()
    CGAL_setup_target_dependencies(${cgal_lib})
  endif()
endforeach()

#
# Define a specific target for basic viewer
#
if (NOT TARGET CGAL::CGAL_Basic_viewer)
  add_library(CGAL::CGAL_Basic_viewer INTERFACE IMPORTED GLOBAL)
    set_target_properties(CGAL::CGAL_Basic_viewer PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "CGAL_USE_BASIC_VIEWER;QT_NO_KEYWORDS"
      INTERFACE_LINK_LIBRARIES CGAL::CGAL_Qt5)
endif()


#warning: the order in this list has to match the enum in Exact_type_selector
set(CGAL_CMAKE_EXACT_NT_BACKEND_OPTIONS GMP_BACKEND GMPXX_BACKEND BOOST_GMP_BACKEND BOOST_BACKEND LEDA_BACKEND MP_FLOAT_BACKEND Default)
set(CGAL_CMAKE_EXACT_NT_BACKEND "Default" CACHE STRING "Setting for advanced users that what to change the default number types used in filtered kernels. Some options might not be working depending on how you configured your build.")
set_property(CACHE CGAL_CMAKE_EXACT_NT_BACKEND PROPERTY STRINGS ${CGAL_CMAKE_EXACT_NT_BACKEND_OPTIONS})

if ( NOT "${CGAL_CMAKE_EXACT_NT_BACKEND}" STREQUAL "Default" )
  list(FIND CGAL_CMAKE_EXACT_NT_BACKEND_OPTIONS ${CGAL_CMAKE_EXACT_NT_BACKEND} DEB_VAL)
  set_property(
      TARGET CGAL
      APPEND PROPERTY
          INTERFACE_COMPILE_DEFINITIONS "CMAKE_OVERRIDDEN_DEFAULT_ENT_BACKEND=${DEB_VAL}"
  ) # do not use set_target_properties to avoid overwritting
endif()
