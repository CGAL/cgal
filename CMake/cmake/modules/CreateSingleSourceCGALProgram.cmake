macro(create_single_source_cgal_program first )

  if(EXISTS ${CMAKE_SOURCE_DIR}/${first})
    set( all ${CMAKE_SOURCE_DIR}/${first} )
    
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_SOURCE_DIR}/${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    add_executable  (${exe_name} ${all})
    
    set_target_properties( ${exe_name} PROPERTIES COMPILE_FLAGS "$(EXTRA_FLAGS) $(TESTSUITE_CXXFLAGS)" )
    set_target_properties( ${exe_name} PROPERTIES LINK_FLAGS    "$(TESTSUITE_LDFLAGS)" )
    
    add_dependencies(${exe_name} CGAL)
    
    if ( CGAL_USE_CGAL_CORE )
      add_dependencies(${exe_name} CGAL_CORE)
    endif()
      
    # Link the executable to CGAL and third-party libraries
    if ( AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} )
    else()
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
    endif()
  
  endif()
    
endmacro()
