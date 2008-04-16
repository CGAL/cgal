# src/Core/CMakeLists.txt
#
# The following variables must be set from the parent script:
# The first two cannot be empty
#
# CGAL_BINARY_DIR
# CGAL_LIB_INSTALL_DIR
#
# CGAL_INCLUDE_DIRS
# GMP_INCLUDE_DIRS
# CGAL_DEFINITIONS
# CGAL_3RD_PARTY_INCLUDE_DIRS
# CGAL_3RD_PARTY_LIBRARIES
# CGAL_3RD_PARTY_DEFINITIONS
# CGAL_LIB_SHARED_LINKER_FLAGS
# CGAL_VERSION
# CGAL_MAJOR_VERSION
# CGAL_MINOR_VERSION
   
assert(CGAL_BINARY_DIR)

add_library(CGAL_CORE "all_files.cpp")

# ${CGAL_CORE_BASENAME} comes defined from FindCGAL_CORE.cmake

if (NOT AUTO_LINK_ENABLED )

  set_target_properties( CGAL_CORE PROPERTIES OUTPUT_NAME ${CGAL_CORE_BASENAME} )
  target_link_libraries( CGAL_CORE ${GMP_LIBRARIES})
  
  set(CGAL_CORE_LIBRARY ${CGAL_CORE_BASENAME} CACHE STRING "The CGAL-core++ library" )
  mark_as_advanced(CGAL_CORE_LIBRARY)
  
else()

  set(CGAL_CORE_LIBRARY)
  
  set(CGAL_CORE_RELEASE_MANGLED_NAME        "${CGAL_CORE_BASENAME}-${TOOLSET}-mt" ) 
  set(CGAL_CORE_DEBUG_MANGLED_NAME          "${CGAL_CORE_BASENAME}-${TOOLSET}-mt-gd" ) 
  set(CGAL_CORE_MINSIZEREL_MANGLED_NAME     "${CGAL_CORE_BASENAME}-${TOOLSET}-mt-o" ) 
  set(CGAL_CORE_RELWITHDEBINFO_MANGLED_NAME "${CGAL_CORE_BASENAME}-${TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL_CORE PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_CORE_DEBUG_MANGLED_NAME}" 
                                              RELEASE_OUTPUT_NAME        "${CGAL_CORE_RELEASE_MANGLED_NAME}" 
                                              MINSIZEREL_OUTPUT_NAME     "${CGAL_CORE_MINSIZEREL_MANGLED_NAME}" 
                                              RELWITHDEBINFO_OUTPUT_NAME "${CGAL_CORE_RELWITHDEBINFO_MANGLED_NAME}" 
                       )
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL_CORE POST_BUILD COMMAND copy /Y "$(TargetPath)" "$(SolutionDir)lib" )
    add_custom_command(TARGET CGAL_CORE POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb" copy /Y "$(TargetDir)$(TargetName).pdb" "$(SolutionDir)lib" )
  endif()
  
endif()

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL_CORE DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()
