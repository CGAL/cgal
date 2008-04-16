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

find_package(ZLIB)
if(ZLIB_FOUND)

  include_directories( ${ZLIB_INCLUDE_DIR} )
  
  add_library(CGAL_IMAGEIO "all_files.cpp")

  target_link_libraries( CGAL_IMAGEIO ${ZLIB_LIBRARY})
  
  set( CGAL_IMAGEIO_BASENAME CGAL-ImageIO )
  
  if (NOT AUTO_LINK_ENABLED )
  
    set_target_properties( CGAL_IMAGEIO PROPERTIES OUTPUT_NAME "${CGAL_IMAGEIO_BASENAME}" )
    
    set(CGAL_IMAGEIO_LIBRARY "${CGAL_IMAGEIO_BASENAME}" CACHE STRING "The CGAL-ImageIO library")
    mark_as_advanced(CGAL_IMAGEIO_LIBRARY)

  else()
  
    set(CGAL_IMAGEIO_LIBRARY)
    
    set(CGAL_IMAGEIO_RELEASE_MANGLED_NAME        "${CGAL_IMAGEIO_BASENAME}-${TOOLSET}-mt" ) 
    set(CGAL_IMAGEIO_DEBUG_MANGLED_NAME          "${CGAL_IMAGEIO_BASENAME}-${TOOLSET}-mt-gd" ) 
    set(CGAL_IMAGEIO_MINSIZEREL_MANGLED_NAME     "${CGAL_IMAGEIO_BASENAME}-${TOOLSET}-mt-o" ) 
    set(CGAL_IMAGEIO_RELWITHDEBINFO_MANGLED_NAME "${CGAL_IMAGEIO_BASENAME}-${TOOLSET}-mt-g" ) 
    
    set_target_properties( CGAL_IMAGEIO PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_IMAGEIO_DEBUG_MANGLED_NAME}" 
                                                   RELEASE_OUTPUT_NAME        "${CGAL_IMAGEIO_RELEASE_MANGLED_NAME}" 
                                                   MINSIZEREL_OUTPUT_NAME     "${CGAL_IMAGEIO_MINSIZEREL_MANGLED_NAME}" 
                                                   RELWITHDEBINFO_OUTPUT_NAME "${CGAL_IMAGEIO_RELWITHDEBINFO_MANGLED_NAME}" 
                        )
    
    if ( HAS_CFG_INTDIR )                      
      add_custom_command(TARGET CGAL_IMAGEIO POST_BUILD COMMAND copy /Y "$(TargetPath)" "$(SolutionDir)lib" )
      add_custom_command(TARGET CGAL_IMAGEIO POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb" copy /Y "$(TargetDir)$(TargetName).pdb" "$(SolutionDir)lib" )
    endif()

  endif()

  if ( CGAL_INSTALL_LIB_DIR )
    install(TARGETS CGAL_IMAGEIO DESTINATION ${CGAL_INSTALL_LIB_DIR} )
  endif()
else()
  message( STATUS "CGAL-ImageIO needs ZLib, cannot be configured.")
endif()  
