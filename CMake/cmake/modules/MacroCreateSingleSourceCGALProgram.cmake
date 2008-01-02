macro(create_single_source_cgal_program sources )

  # Create an executable for each .cpp in folder
  foreach(CPP ${sources})

    message( STATUS "CGAL program: ${CPP}"  )
    
    # Add executable that is built from the source file ${CPP}.
    # The executable is named after ${CPP}'s basename. The extensions are automatically found.
    get_filename_component(exe_name ${CPP} NAME_WE)
    
    add_executable (${exe_name} ${CPP})
    
    # Link the executable to CGAL and third-party libraries
    if ( NOT AUTO_LINK_ENABLED )    
      target_link_libraries(${exe_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    endif()

    # Add " make test" rule for executable
    # TODO: get test parameters from ${exe_name}.cmd
    add_test(${exe_name} ${exe_name})
    
  endforeach()
  
endmacro()
