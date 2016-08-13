function(create_single_source_cgal_program firstfile )

  if(NOT IS_ABSOLUTE "${firstfile}")
    set(firstfile "${CMAKE_CURRENT_SOURCE_DIR}/${firstfile}")
  endif()

  get_filename_component(exe_name ${firstfile} NAME_WE)

  if(EXISTS "${firstfile}")

    set( all "${firstfile}" )

    # remaining files
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_CURRENT_SOURCE_DIR}/${i} )
    endforeach()


    add_executable(${exe_name} ${all})

    if(BUILD_TESTING)
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cin")
	add_test(NAME ${exe_name}
          COMMAND ${CMAKE_COMMAND}
	    -DCMD:STRING=$<TARGET_FILE:${exe_name}>
	    -DCIN:STRING=${exe_name}.cin
	    -P ${CGAL_MODULES_DIR}/run_test_with_cin.cmake
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	message(STATUS "add test: ${exe_name} < ${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cin")
      else()
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cmd")
          file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cmd" ARGS)
	endif()
	message(STATUS "add test: ${exe_name} ${ARGS}")
	add_test(NAME ${exe_name}
          COMMAND ${exe_name} ${ARGS}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      endif()
      set_property(TEST "${exe_name}"
        APPEND PROPERTY LABELS "${PROJECT_NAME}")
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
