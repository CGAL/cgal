# This module install a hook (with `variable_watch()`) on CMAKE_CURRENT_LIST_DIR
# 
# That uses the non-documented fact that CMAKE_CURRENT_LIST_DIR is cleared
# by CMake at the end of the configuration process. So, if the value of
# that variable gets empty, that means that CMake has reached the end of
# the configuration.
#
# See https://stackoverflow.com/a/43300621/1728537 for the starting point.

get_property(PROPERTY_CGAL_run_at_the_end_of_configuration_INCLUDED
  GLOBAL PROPERTY CGAL_run_at_the_end_of_configuration_INCLUDED)
if(PROPERTY_CGAL_run_at_the_end_of_configuration_INCLUDED)
  return()
endif()

function(CGAL_run_at_the_end_of_configuration variable access value current_list_file stack)
  if(NOT access STREQUAL "MODIFIED_ACCESS" OR value)
    # Only do something at the end of the CMake process, when the value of
    # variable CMAKE_CURRENT_LIST_DIR is changed to the empty string.
    return()
  endif()
  # Warn when CMAKE_BUILD_TYPE is empty or Debug
  if(DEFINED CMAKE_BUILD_TYPE AND ( NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug") )
    set(keyword WARNING)
    set(type warning)
    if(RUNNING_CGAL_AUTO_TEST)
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

variable_watch("CMAKE_CURRENT_LIST_DIR" CGAL_run_at_the_end_of_configuration)

set_property(GLOBAL PROPERTY CGAL_run_at_the_end_of_configuration_INCLUDED TRUE)
