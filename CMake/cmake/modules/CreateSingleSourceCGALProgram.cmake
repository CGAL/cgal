macro(create_single_source_cgal_program first )

  if(EXISTS ${first})
    set( all ${first} )
    
    foreach( i ${ARGN} )
      set( all ${all} ${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    add_executable  (${exe_name} ${all})
    
    #add_dependencies(${exe_name} CGAL CGAL_CORE)
      
    # Link the executable to CGAL and third-party libraries
    if ( NOT AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    endif()
  
    # TODO: get test parameters from ${exe_name}.cmd
    #set ( test_args )
    # 
    #if ( EXISTS ${exe_name}.cmd )
    #  file( READ ${exe_name}.cmd ${test_args} )
    #  message( STATUS "Command line arguments for ${exe_name}: ${test_args} )  
    #endif()
    # 
    #add_test( ${exe_name} ${exe_name} ${ARGS} )
  endif(EXISTS ${first})
    
endmacro()
