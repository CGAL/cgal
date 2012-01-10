build_cgal_library(CGAL CGAL "") 

foreach (lib ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES})

  set (vlib "${CGAL_EXT_LIB_${lib}_PREFIX}")

  if (CGAL_USE_${vlib})

    message( STATUS "External library ${vlib} has been preconfigured")
    use_lib( ${vlib} "###${${vlib}_USE_FILE}")

  endif()

endforeach()

include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

set( CGAL_BASENAME CGAL )

set_target_properties(CGAL PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")

if (NOT CGAL_AUTO_LINK_ENABLED )

  set_target_properties( CGAL PROPERTIES OUTPUT_NAME "${CGAL_BASENAME}" )
  
  get_target_property( CGAL_LIBRARY_FULL_PATH CGAL LOCATION  )
  
  set(CGAL_LIBRARY "${CGAL_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL library" FORCE)
  
  hide_variable(CGAL_LIBRARY)
  
else()
  
  set(CGAL_Release_MANGLED_NAME        "${CGAL_BASENAME}-${CGAL_TOOLSET}-mt" ) 
  set(CGAL_Debug_MANGLED_NAME          "${CGAL_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
  set(CGAL_MinSizeRel_MANGLED_NAME     "${CGAL_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
  set(CGAL_RelWithDebInfo_MANGLED_NAME "${CGAL_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_Debug_MANGLED_NAME}" 
                                         RELEASE_OUTPUT_NAME        "${CGAL_Release_MANGLED_NAME}" 
                                         MINSIZEREL_OUTPUT_NAME     "${CGAL_MinSizeRel_MANGLED_NAME}" 
                                         RELWITHDEBINFO_OUTPUT_NAME "${CGAL_RelWithDebInfo_MANGLED_NAME}" 
                       )
  
  
  set(CGAL_LIBRARY "${CGAL_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL library" FORCE)
  set(CGAL_LIBRARY "${CGAL_${CMAKE_BUILD_TYPE}_MANGLED_NAME}")
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL POST_BUILD COMMAND copy /Y \"$(TargetPath)\" \"$(SolutionDir)lib\" )
    add_custom_command(TARGET CGAL POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(SolutionDir)lib\" )
  endif()
  
endif()

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

