# src/CGAL/CMakeLists.txt
#
# The following variables must be set from the parent script:
# The first two cannot be empty
#
# CGAL_BINARY_DIR
# CGAL_LIB_INSTALL_DIR
#
# CGAL_INCLUDE_DIRS
# CGAL_DEFINITIONS
# CGAL_3RD_PARTY_INCLUDE_DIRS
# CGAL_3RD_PARTY_LIBRARIES
# CGAL_3RD_PARTY_DEFINITIONS
# CGAL_LIB_SHARED_LINKER_FLAGS
# CGAL_VERSION
# CGAL_MAJOR_VERSION
# CGAL_MINOR_VERSION

   
assert(CGAL_BINARY_DIR)

add_library (CGAL "all_files.cpp" )

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

if (CGAL_USE_CGAL_CORE)
  add_dependencies(CGAL CGAL_CORE)
endif()

set( CGAL_BASENAME CGAL )

set_target_properties(CGAL PROPERTIES VERSION "${CGAL_VERSION}" SOVERSION "${CGAL_SOVERSION}")

if (NOT AUTO_LINK_ENABLED )

  set_target_properties( CGAL PROPERTIES OUTPUT_NAME "${CGAL_BASENAME}" )
  
  get_target_property( CGAL_LIBRARY_FULL_PATH CGAL LOCATION  )
  
  set(CGAL_LIBRARY "${CGAL_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL library" )
  
  hide_variable(CGAL_LIBRARY)
  
else()

  set(CGAL_LIBRARY)
  
  set(CGAL_RELEASE_MANGLED_NAME        "${CGAL_BASENAME}-${TOOLSET}-mt" ) 
  set(CGAL_DEBUG_MANGLED_NAME          "${CGAL_BASENAME}-${TOOLSET}-mt-gd" ) 
  set(CGAL_MINSIZEREL_MANGLED_NAME     "${CGAL_BASENAME}-${TOOLSET}-mt-o" ) 
  set(CGAL_RELWITHDEBINFO_MANGLED_NAME "${CGAL_BASENAME}-${TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_DEBUG_MANGLED_NAME}" 
                                         RELEASE_OUTPUT_NAME        "${CGAL_RELEASE_MANGLED_NAME}" 
                                         MINSIZEREL_OUTPUT_NAME     "${CGAL_MINSIZEREL_MANGLED_NAME}" 
                                         RELWITHDEBINFO_OUTPUT_NAME "${CGAL_RELWITHDEBINFO_MANGLED_NAME}" 
                       )
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL POST_BUILD COMMAND copy /Y "$(TargetPath)" "$(SolutionDir)lib" )
    add_custom_command(TARGET CGAL POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb" copy /Y "$(TargetDir)$(TargetName).pdb" "$(SolutionDir)lib" )
  endif()
  
endif()

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

