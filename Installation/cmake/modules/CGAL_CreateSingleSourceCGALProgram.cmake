function(create_single_source_cgal_program first )

  set(file "${first}")

  if(NOT IS_ABSOLUTE "${file}")
    set(file "${CMAKE_CURRENT_SOURCE_DIR}/${first}")
  endif()

  if(EXISTS "${file}")
  
    set( all "${file}" )
    
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_CURRENT_SOURCE_DIR}/${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    add_executable(${exe_name} ${all})
    
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${exe_name} )
    
    # Link the executable to CGAL and third-party libraries
    if ( CGAL_AUTO_LINK_ENABLED )    

      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} )

    else()

      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

    endif()
  
  endif()
    
endfunction()
