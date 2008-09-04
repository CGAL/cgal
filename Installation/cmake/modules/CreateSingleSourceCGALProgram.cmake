macro(create_single_source_cgal_program prefix first )

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${first})
  
    set( all ${CMAKE_CURRENT_SOURCE_DIR}/${first} )
    
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_CURRENT_SOURCE_DIR}/${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    set( target_name "${prefix}_${exe_name}" )
    
    add_executable(${target_name} ${all})
    
    set( CGAL_EXECUTABLE_TARGETS )
    set( CGAL_EXECUTABLE_TARGETS "${CGAL_EXECUTABLE_TARGETS}" "${target_name}" CACHE INTERNAL "" FORCE )
    
    # Link the executable to CGAL and third-party libraries
    if ( AUTO_LINK_ENABLED )    
      target_link_libraries(${target_name} ${CGAL_3RD_PARTY_LIBRARIES} )
    else()
      target_link_libraries(${target_name} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_LIBRARIES} )
    endif()
  
  endif()
    
endmacro()
