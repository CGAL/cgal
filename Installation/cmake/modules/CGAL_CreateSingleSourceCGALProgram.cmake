function(create_single_source_cgal_program firstfile )

  if(NOT IS_ABSOLUTE "${firstfile}")
    set(firstfile "${CGAL_CURRENT_SOURCE_DIR}/${firstfile}")
  endif()

  get_filename_component(exe_name ${firstfile} NAME_WE)

  if(EXISTS "${firstfile}")

    set( all "${firstfile}" )

    # remaining files
    foreach( i ${ARGN} )
      set( all ${all} ${CGAL_CURRENT_SOURCE_DIR}/${i} )
    endforeach()


    add_executable(${exe_name} ${all})

    if(BUILD_TESTING)
      set(cin_file "${CGAL_CURRENT_SOURCE_DIR}/${exe_name}.cin")
      if(EXISTS ${cin_file})
	add_test(NAME ${exe_name}
          COMMAND ${CMAKE_COMMAND}
	    -DCMD:STRING=$<TARGET_FILE:${exe_name}>
	    -DCIN:STRING=${cin_file}
	    -P "${CGAL_MODULES_DIR}/run_test_with_cin.cmake")
#	message(STATUS "add test: ${exe_name} < ${cin_file}")
      else()
	# TODO: deal with shell globbing; if the `cmd` file contains
	# a `*`, then interprete the command using bash.
	set(cmd_file "${CGAL_CURRENT_SOURCE_DIR}/${exe_name}.cmd")
	if(EXISTS ${cmd_file})
          file(STRINGS "${cmd_file}" CMD_LINES)
	  set(ARGS)
#	  message(STATUS "DEBUG test ${exe_name}")
	  foreach(CMD_LINE ${CMD_LINES})
#	    message(STATUS "  command line: ${CMD_LINE}")
	    separate_arguments(CMD_LINE_ARGS UNIX_COMMAND ${CMD_LINE})
#	    message(STATUS "  args: ${CMD_LINE_ARGS}")
	    list(APPEND ARGS ${CMD_LINE_ARGS})
	  endforeach()
	endif()
#	message(STATUS "add test: ${exe_name} ${ARGS}")
	add_test(NAME ${exe_name} COMMAND ${exe_name} ${ARGS})
      endif()
      set_property(TEST "${exe_name}"
        APPEND PROPERTY LABELS "${PROJECT_NAME}")
#      message(STATUS "  workding dir: ${CGAL_CURRENT_SOURCE_DIR}")
      set_property(TEST "${exe_name}"
        PROPERTY WORKING_DIRECTORY ${CGAL_CURRENT_SOURCE_DIR})
    endif(BUILD_TESTING)

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${exe_name} )

    # Link the executable to CGAL and third-party libraries
    if ( CGAL_AUTO_LINK_ENABLED )

      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} )

    else()

      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

    endif()
  else()
    message(AUTHOR_WARNING "The executable ${exe_name} will not be created because the source file ${firstfile} does not exist.")
  endif()

endfunction()
