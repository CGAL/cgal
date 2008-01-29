macro(create_single_source_cgal_program first )

  set( all ${first} )
  
  foreach( i ${ARGN} )
    set( all ${all} ${i} ) 
  endforeach()
  
  get_filename_component(exe_name ${first} NAME_WE)
  
  add_executable  (${exe_name} ${all})
  add_dependencies(${exe_name} CGAL CGAL_CORE)
    
  # Link the executable to CGAL and third-party libraries
  if ( NOT AUTO_LINK_ENABLED )    
    target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  endif()

  # Add " make test" rule for executable
  # TODO: get test parameters from ${exe_name}.cmd
  # add_test(${exe_name} ${exe_name})
  
endmacro()
