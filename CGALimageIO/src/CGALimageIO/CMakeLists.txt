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
find_package(OpenGL)

option(WITH_VTK "Add VTK support to the CGAL ImageIO library." OFF)

if(WITH_VTK)
  find_package(VTK QUIET)
endif(WITH_VTK)

if(ZLIB_FOUND AND OPENGL_FOUND)

  get_dependency_version(ZLIB)
  get_dependency_version(OPENGL OpenGL)
   
  message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
  message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
  message( STATUS "OpenGL definitions: ${OPENGL_DEFINITIONS}" )
  message( STATUS "ZLib   include:     ${ZLIB_INCLUDE_DIR}" )
  message( STATUS "ZLib   libraries:   ${ZLIB_LIBRARY}" )
  
  set( CGAL_ImageIO_BASENAME CGAL-ImageIO )
  
  cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
  cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES    ${ZLIB_LIBRARY} ${OPENGL_LIBRARIES} )
  cache_set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS  "" )
  
  if(WITH_VTK)
  
    if( VTK_FOUND )
    
      message(STATUS "VTK-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION} found. Build VTK support in ${CGAL_ImageIO_BASENAME}.")
      include( ${VTK_USE_FILE} )
      cache_set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS  ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS} -DCGAL_USE_VTK )
      cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} ${VTK_INCLUDE_DIR} )
      cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES    ${CGAL_ImageIO_3RD_PARTY_LIBRARIES} vtkImaging vtkIO )
      
    else()
      message(STATUS "VTK not found.")
    endif()
    
  endif()
  
  include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} )
  
  add_library(CGAL_ImageIO all_files.cpp)
  
  add_dependencies( CGAL_ImageIO CGAL )
  
  add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS} )
  
  set_target_properties(CGAL_ImageIO PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
  
  target_link_libraries( CGAL_ImageIO ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES} )
  
  if (NOT AUTO_LINK_ENABLED )
  
    set_target_properties( CGAL_ImageIO PROPERTIES OUTPUT_NAME "${CGAL_ImageIO_BASENAME}" )
    
    get_target_property( CGAL_ImageIO_LIBRARY_FULL_PATH CGAL_ImageIO LOCATION  )
    
    set(CGAL_ImageIO_LIBRARY "${CGAL_ImageIO_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL-ImageIO library")
    
    mark_as_advanced(CGAL_ImageIO_LIBRARY)

    target_link_libraries( CGAL_ImageIO ${CGAL_LIBRARY} )

  else()
  
    set(CGAL_ImageIO_LIBRARY)
    
    set(CGAL_ImageIO_RELEASE_MANGLED_NAME        "${CGAL_ImageIO_BASENAME}-${TOOLSET}-mt" ) 
    set(CGAL_ImageIO_DEBUG_MANGLED_NAME          "${CGAL_ImageIO_BASENAME}-${TOOLSET}-mt-gd" ) 
    set(CGAL_ImageIO_MINSIZEREL_MANGLED_NAME     "${CGAL_ImageIO_BASENAME}-${TOOLSET}-mt-o" ) 
    set(CGAL_ImageIO_RELWITHDEBINFO_MANGLED_NAME "${CGAL_ImageIO_BASENAME}-${TOOLSET}-mt-g" ) 
    
    set_target_properties( CGAL_ImageIO PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_ImageIO_DEBUG_MANGLED_NAME}" 
                                                   RELEASE_OUTPUT_NAME        "${CGAL_ImageIO_RELEASE_MANGLED_NAME}" 
                                                   MINSIZEREL_OUTPUT_NAME     "${CGAL_ImageIO_MINSIZEREL_MANGLED_NAME}" 
                                                   RELWITHDEBINFO_OUTPUT_NAME "${CGAL_ImageIO_RELWITHDEBINFO_MANGLED_NAME}" 
                        )
    
    if ( HAS_CFG_INTDIR )                      
      add_custom_command(TARGET CGAL_ImageIO POST_BUILD COMMAND copy /Y "$(TargetPath)" "$(SolutionDir)lib" )
      add_custom_command(TARGET CGAL_ImageIO POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb" copy /Y "$(TargetDir)$(TargetName).pdb" "$(SolutionDir)lib" )
    endif()

  endif()

  if ( CGAL_INSTALL_LIB_DIR )
    install(TARGETS CGAL_ImageIO DESTINATION ${CGAL_INSTALL_LIB_DIR} )
  endif()
else()
  message( STATUS "CGAL-ImageIO needs ZLib and OpenGL, cannot be configured.")
endif() 
