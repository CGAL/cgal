macro(create_single_source_cgal_qt4_program first )

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${first})
    set( all ${CMAKE_CURRENT_SOURCE_DIR}/${first} )
  
    foreach( i ${ARGN} )
      set( all ${all} ${CMAKE_CURRENT_SOURCE_DIR}/${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    QT4_AUTOMOC( ${all} )
    
    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
    
    add_executable  (${exe_name} ${all})
    
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${exe_name} )
    
    # Link the executable to CGAL and third-party libraries
    if ( AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
    else()
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
    endif()

  endif(EXISTS ${first})
  
endmacro()
