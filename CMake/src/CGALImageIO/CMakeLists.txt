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
assert(CGAL_LIB_INSTALL_DIR)

find_package(ZLIB)
if(ZLIB_FOUND)

  include_directories( ${ZLIB_INCLUDE_DIR} )
  
  add_library(CGAL_IMAGE_IO "all_files.cpp")

  target_link_libraries( CGAL_IMAGE_IO ${ZLIB_LIBRARY})
  
  if (NOT AUTO_LINK_ENABLED )
  
    set_target_properties( CGAL_IMAGE_IO PROPERTIES OUTPUT_NAME "cgal-image-io" )
    
  else()
  
    set(CGAL_IMAGE_IO_RELEASE_MANGLED_NAME        "cgal-image-io${TOOLSET}-mt" ) 
    set(CGAL_IMAGE_IO_DEBUG_MANGLED_NAME          "cgal-image-io${TOOLSET}-mt-gd" ) 
    set(CGAL_IMAGE_IO_MINSIZEREL_MANGLED_NAME     "cgal-image-io${TOOLSET}-mt-o" ) 
    set(CGAL_IMAGE_IO_RELWITHDEBINFO_MANGLED_NAME "cgal-image-io${TOOLSET}-mt-g" ) 
    
    set_target_properties( CGAL_IMAGE_IO PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_IMAGE_IO_DEBUG_MANGLED_NAME}" 
                                                    RELEASE_OUTPUT_NAME        "${CGAL_IMAGE_IO_RELEASE_MANGLED_NAME}" 
                                                    MINSIZEREL_OUTPUT_NAME     "${CGAL_IMAGE_IO_MINSIZEREL_MANGLED_NAME}" 
                                                    RELWITHDEBINFO_OUTPUT_NAME "${CGAL_IMAGE_IO_RELWITHDEBINFO_MANGLED_NAME}" 
                        )
    
    add_custom_command(TARGET CGAL_IMAGE_IO POST_BUILD COMMAND move /Y $(TargetDir)$(TargetName)$(TargetExt) ../../lib )
    add_custom_command(TARGET CGAL_IMAGE_IO POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb." move /Y $(TargetDir)$(TargetName).pdb ../../lib)
  endif()

  # WARNING: Use only relative paths; full paths break CPack!
  install(TARGETS CGAL_IMAGE_IO DESTINATION ${CGAL_LIB_INSTALL_DIR} )
else()
  message( WARNING "cgal-image-io needs ZLib, cannot be configured.")
endif()  
