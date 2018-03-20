function(CGAL_run_at_the_end_of_configuration variable access value current_list_file stack)
  if(value)
    # Only do something at the end of the CMake process
    return()
  endif()
  # Warn when CMAKE_BUILD_TYPE is empty or Debug
  if(DEFINED CMAKE_BUILD_TYPE AND ( NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug") )
    set(keyword WARNING)
    if(RUNNING_CGAL_AUTO_TEST)
      # No warning in the CMake test suite, but a status message
      set(keyword )
    endif()
    if(NOT CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE)
      message(${keyword} "\
=======================================================================\n\
CGAL performance notice:\n\
The variable CMAKE_BUILD_TYPE is set to \"${CMAKE_BUILD_TYPE}\".  That \
value should be used only for debugging. For performance reasons, you \
should set CMAKE_BUILD_TYPE to \"Release\".\n\
Set CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE to TRUE if you want to \
disable this warning.\n\
=======================================================================\
")
    endif()
  endif()
endfunction()

variable_watch("CMAKE_CURRENT_LIST_DIR" CGAL_run_at_the_end_of_configuration)
