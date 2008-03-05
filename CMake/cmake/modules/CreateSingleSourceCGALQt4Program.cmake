macro(create_single_source_cgal_qt4_program first )

  if(EXISTS ${first})
    set( all ${first} )
  
    foreach( i ${ARGN} )
      set( all ${all} ${i} ) 
    endforeach()
    
    get_filename_component(exe_name ${first} NAME_WE)
    
    QT4_AUTOMOC( ${all} )
    
    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
    
    add_executable  (${exe_name} ${all})
    add_dependencies(${exe_name} CGAL CGAL_CORE)
    
    set_target_properties( ${exe_name} PROPERTIES COMPILE_FLAGS "$(EXTRA_FLAGS) $(TESTSUITE_CXXFLAGS)" )
    set_target_properties( ${exe_name} PROPERTIES LINK_FLAGS    "$(TESTSUITE_LDFLAGS)" )
    
    # Link the executable to CGAL and third-party libraries
    if ( AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
    else()
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_QT_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES})
    endif()

  endif(EXISTS ${first})
  
endmacro()
