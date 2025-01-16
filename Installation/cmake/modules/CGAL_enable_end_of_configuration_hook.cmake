# This module install a hook (with `variable_watch()`) on CMAKE_CURRENT_LIST_DIR
#
# That uses the non-documented fact that CMAKE_CURRENT_LIST_DIR is cleared
# by CMake at the end of the configuration process. So, if the value of
# that variable gets empty, that means that CMake has reached the end of
# the configuration.
#
# See https://stackoverflow.com/a/43300621/1728537 for the starting point.
cmake_minimum_required(VERSION 3.12...3.29)

if(CGAL_SKIP_CMAKE_HOOKS)
  return()
endif()
get_property(PROPERTY_CGAL_run_at_the_end_of_configuration_INCLUDED
  GLOBAL PROPERTY CGAL_run_at_the_end_of_configuration_INCLUDED)
if(PROPERTY_CGAL_run_at_the_end_of_configuration_INCLUDED)
  CGAL_install_hooks()
  return()
endif()

function(CGAL_hook_check_targets)
  if(NOT CGAL_CHECK_UNREFERENCES_TARGETS)
    return()
  endif()
  get_property(_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
  set(_list_of_deps)
  set(_special_targets  demos examples tests ALL_CGAL_TARGETS CGAL_Qt6_moc_and_resources uninstall install_FindCGAL)
  foreach(t ${_special_targets})
    if(NOT TARGET ${t})
      continue()
    endif()
    get_property(_deps TARGET ${t} PROPERTY MANUALLY_ADDED_DEPENDENCIES)
    # message("  deps of ${t}: ${_deps}")
    list(APPEND _list_of_deps ${_deps})
  endforeach()
  list(APPEND _list_of_deps ${CGAL_EXECUTABLE_TARGETS})
  # message(STATUS "  all deps: ${_list_of_deps}")
  foreach(target ${_targets})
    # message(STATUS "  new target: ${target}")
    if(${target} IN_LIST _special_targets)
      continue()
    endif()
    if(NOT ${target} IN_LIST _list_of_deps)
      message(AUTHOR_WARNING "  orphan target: ${target}")
    endif()
  endforeach()
endfunction()

function(CGAL_hook_check_unused_cpp_files)
  if(NOT CGAL_CHECK_UNUSED_CPP_FILES)
    return()
  endif()
  get_property(_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
  file(GLOB _cppfiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  if(NOT _targets OR NOT _cppfiles)
    return()
  endif()
  set(_sources)
  foreach(_target ${_targets})
    get_property(_target_type TARGET ${_target} PROPERTY TYPE)
    if(_target_type STREQUAL INTERFACE_LIBRARY)
      continue()
    endif()
    get_property(_target_sources TARGET ${_target} PROPERTY SOURCES)
    list(APPEND _sources ${_target_sources})
  endforeach()
  if(_sources)
    list(REMOVE_ITEM _cppfiles ${_sources})
  endif()
  if(_cppfiles)
    set(_warning "In ${CMAKE_CURRENT_SOURCE_DIR}, the following files are unused:")
    foreach(_cppfile ${_cppfiles})
      set(_warning "${_warning}
  ${_cppfile}")
    endforeach()
    set(_warning  "${_warning}
  ")
    message(AUTHOR_WARNING "${_warning}")
  endif()
endfunction()

function(CGAL_hook_check_CMAKE_BUILD_TYPE)
  # Warn when CMAKE_BUILD_TYPE is empty or Debug
  if(DEFINED CMAKE_BUILD_TYPE AND ( NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug") )
    set(keyword WARNING)
    set(type warning)
    if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
      # No warning in the CMake test suite, but a status message
      set(keyword)
      set(type notice)
    endif()
    if(NOT CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE)
      message(${keyword} "\
=======================================================================\n\
CGAL performance notice:\n\
The variable CMAKE_BUILD_TYPE is set to \"${CMAKE_BUILD_TYPE}\". For \
performance reasons, you should set CMAKE_BUILD_TYPE to \"Release\".\n\
Set CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE to TRUE if you want to \
disable this ${type}.\n\
=======================================================================\
")
    endif()
  endif()
endfunction()

function(CGAL_hook_fix_ctest_depending_on_Qt6)
  get_property(_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target ${_targets})
    if(NOT TEST "compilation of  ${_target}")
      continue()
    endif()
    get_property(_target_links TARGET ${_target} PROPERTY LINK_LIBRARIES)
    if("CGAL_Qt6" IN_LIST _target_links OR "CGAL::CGAL_Qt6" IN_LIST _target_links OR "CGAL::CGAL_Basic_viewer" IN_LIST _target_links)
      set_property(TEST "compilation of  ${_target}" APPEND PROPERTY FIXTURES_REQUIRED CGAL_Qt6_moc_and_resources_Fixture)
    endif()
    endforeach()
endfunction()

function(CGAL_hook_fix_ctest_dependencies)
  get_property(_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
  # message(STATUS "  ${CMAKE_CURRENT_SOURCE_DIR} targets: ${_targets}")
  foreach(_target ${_targets})
    if(NOT TEST "compilation of  ${_target}")
      continue()
    endif()
    get_property(linked_libraries TARGET "${_target}" PROPERTY LINK_LIBRARIES)
    get_property(deps TARGET "${_target}" PROPERTY MANUALLY_ADDED_DEPENDENCIES)
    list(APPEND deps ${linked_libraries})
    # message("  ${_target} depends on: ${deps}")
    foreach(dep ${deps})
      if(TARGET ${dep} AND NOT dep MATCHES "^CGAL::|^CGAL_|^CGAL$" )
        get_target_property(imported ${dep} IMPORTED)
        if(NOT imported)
          # message(STATUS "The target ${dep} is a dependency of ${_target}.")
          set_property(TEST "compilation of  ${_target}"
            APPEND PROPERTY DEPENDS "compilation of  ${dep}")
        endif()
      endif()
    endforeach()
  endforeach()
endfunction()

function(CGAL_hooks_at_end_of_all_directories)
  CGAL_hook_check_targets()
  CGAL_hook_check_unused_cpp_files()
  if(CGAL_ENABLE_TESTING)
    CGAL_hook_fix_ctest_depending_on_Qt6()
    CGAL_hook_fix_ctest_dependencies()
  endif()
endfunction()

function(CGAL_run_at_the_end_of_configuration variable access value current_list_file stack)
  if(NOT access STREQUAL "MODIFIED_ACCESS")
    return()
  endif()
  if(current_list_file STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
        AND NOT current_list_file STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
        OR NOT stack STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
        AND stack STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
        OR stack MATCHES doc/CMakeLists.txt
        OR stack MATCHES demo/Polyhedron/)
      return()
  endif()


  if(NOT current_list_file STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
     AND stack STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    CGAL_hooks_at_end_of_all_directories()
  endif()

  if(value)
    # Only do the following at the end of the CMake process, when the
    # value of variable CMAKE_CURRENT_LIST_DIR is changed to the empty
    # string.
    return()
  endif()

  CGAL_hook_check_CMAKE_BUILD_TYPE()

endfunction()

function(CGAL_install_hooks)
  get_property(PROPERTY_CGAL_hooks_installed DIRECTORY PROPERTY CGAL_hooks_installed)
  if(PROPERTY_CGAL_hooks_installed)
    return()
  endif()
  # message(STATUS "Installing hooks in ${CMAKE_CURRENT_SOURCE_DIR}")
  if(CMAKE_VERSION VERSION_LESS 3.19)
    variable_watch("CMAKE_CURRENT_LIST_DIR" CGAL_run_at_the_end_of_configuration)
  else()
    cmake_language(DEFER CALL CGAL_hooks_at_end_of_all_directories)
    if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
      cmake_language(DEFER CALL CGAL_hook_check_CMAKE_BUILD_TYPE)
    endif()
  endif()
  set_property(DIRECTORY PROPERTY CGAL_hooks_installed TRUE)
endfunction()

CGAL_install_hooks()

set_property(GLOBAL PROPERTY CGAL_run_at_the_end_of_configuration_INCLUDED TRUE)
