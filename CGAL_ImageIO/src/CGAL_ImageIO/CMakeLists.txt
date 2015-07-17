message("Configuring libCGAL_ImageIO")

find_package( OpenGL )
find_package( ZLIB )

#option(WITH_VTK "Add VTK support to libCGAL_ImageIO." OFF)

if(WITH_VTK)
  find_package(VTK COMPONENTS vtkImagingCore vtkIOImage vtkFiltersImaging NO_MODULE REQUIRED)
endif(WITH_VTK)

if(OPENGL_FOUND)
  message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
  message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
  
  cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR} )
  cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES    ${OPENGL_LIBRARIES}   )

  if(ZLIB_FOUND)
    get_dependency_version(ZLIB)
    cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR} )
    cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES    ${CGAL_ImageIO_3RD_PARTY_LIBRARIES} ${ZLIB_LIBRARY} )
    cache_set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS  "-DCGAL_USE_ZLIB")
    cache_set(CGAL_ImageIO_USE_ZLIB "ON")
  endif(ZLIB_FOUND)

  set( CGAL_ImageIO_BASENAME CGAL_ImageIO)
  
  if(WITH_VTK)
    if( VTK_FOUND )
      
      message(STATUS "VTK-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION} found. Build VTK support in ${CGAL_ImageIO_BASENAME}.")
      include( ${VTK_USE_FILE} )
      cache_set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS  ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS} -DCGAL_USE_VTK)
      cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS} )
      cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES_DIRS ${CGAL_ImageIO_3RD_PARTY_LIBRARIES_DIRS} ${vtkImagingCore_RUNTIME_LIBRARY_DIRS} ${vtkImagingCore_LIBRARY_DIRS} ${vtkIOImage_RUNTIME_LIBRARY_DIRS} ${vtkIOImage_LIBRARY_DIRS} ${vtkFiltersImaging_RUNTIME_LIBRARY_DIRS} ${vtkFiltersImaging_LIBRARY_DIRS})
      cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES    ${CGAL_ImageIO_3RD_PARTY_LIBRARIES}   ${VTK_LIBRARIES}  )
    else()
      message(STATUS "VTK not found.")
    endif()
  endif()
  
  if(COMMAND add_config_flag)
    set( CGAL_HAS_IMAGEIO TRUE ) 
    add_config_flag( CGAL_HAS_IMAGEIO )
  endif()

  use_essential_libs()

  include_directories( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} )
  
  link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )

  collect_cgal_library( CGAL_ImageIO "")
  
  add_dependencies( CGAL_ImageIO CGAL )
  
  add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS} )
  
  # CGAL_ImageIO only depends on CGAL in DEBUG, but we still link it
  # in both build types.
  target_link_libraries( CGAL_ImageIO CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES} )
  
message("libCGAL_ImageIO is configured")
    
else(OPENGL_FOUND)
  message( STATUS "libCGAL_ImageIO needs OpenGL, cannot be configured.")
endif(OPENGL_FOUND)
  
if(NOT ZLIB_FOUND)
  message( STATUS "NOTICE: libCGAL_ImageIO needs ZLib to read compressed files. That feature will not be activated.")
endif(NOT ZLIB_FOUND) 

