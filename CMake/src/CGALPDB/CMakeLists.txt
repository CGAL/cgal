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

add_library(CGAL_PDB "all_files.cpp")

set( CGAL_PDB_BASENAME CGAL-PDB )

if (NOT AUTO_LINK_ENABLED )

  set_target_properties( CGAL_PDB PROPERTIES OUTPUT_NAME "${CGAL_PDB_BASENAME}" )
  
  set(CGAL_PDB_LIBRARY "${CGAL_PDB_BASENAME}" CACHE STRING "The CGAL-PDB library")
  mark_as_advanced(CGAL_PDB_LIBRARY)
  
else()

  set(CGAL_PDB_LIBRARY)
  
  set(CGAL_PDB_RELEASE_MANGLED_NAME        "${CGAL_PDB_BASENAME}-${TOOLSET}-mt" ) 
  set(CGAL_PDB_DEBUG_MANGLED_NAME          "${CGAL_PDB_BASENAME}-${TOOLSET}-mt-gd" ) 
  set(CGAL_PDB_MINSIZEREL_MANGLED_NAME     "${CGAL_PDB_BASENAME}-${TOOLSET}-mt-o" ) 
  set(CGAL_PDB_RELWITHDEBINFO_MANGLED_NAME "${CGAL_PDB_BASENAME}-${TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL_PDB PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_PDB_DEBUG_MANGLED_NAME}" 
                                              RELEASE_OUTPUT_NAME        "${CGAL_PDB_RELEASE_MANGLED_NAME}" 
                                              MINSIZEREL_OUTPUT_NAME     "${CGAL_PDB_MINSIZEREL_MANGLED_NAME}" 
                                              RELWITHDEBINFO_OUTPUT_NAME "${CGAL_PDB_RELWITHDEBINFO_MANGLED_NAME}" 
                       )
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL_PDB POST_BUILD COMMAND copy /Y $(TargetDir)$(TargetName)$(TargetExt) ../../lib )
    add_custom_command(TARGET CGAL_PDB POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb." copy /Y $(TargetDir)$(TargetName).pdb ../../lib)
  endif()
  
endif()

# WARNING: Use only relative paths; full paths break CPack!
install(TARGETS CGAL_PDB DESTINATION ${CGAL_INSTALL_LIB_DIR} )
