message("Configuring libCGAL_ImageIO")

if(NOT CGAL_HEADER_ONLY)
  find_package( OpenGL )
  find_package( ZLIB )

  if(OPENGL_FOUND)
    message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
    message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )

    cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS "${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS}" ${OPENGL_INCLUDE_DIR})
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

    include_directories( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS} .)

    link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
  endif(OPENGL_FOUND)

else (NOT CGAL_HEADER_ONLY)
  cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(CGAL_HEADER_ONLY OR OPENGL_FOUND)
  collect_cgal_library( CGAL_ImageIO "")
endif()

if(NOT CGAL_HEADER_ONLY AND OPENGL_FOUND)
  add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS} )

  add_dependencies( CGAL_ImageIO CGAL )

  # CGAL_ImageIO only depends on CGAL in DEBUG, but we still link it
  # in both build types.
  target_link_libraries( CGAL_ImageIO CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES} )
endif()

if(NOT CGAL_HEADER_ONLY AND NOT OPENGL_FOUND)
  message( STATUS "libCGAL_ImageIO needs OpenGL, cannot be configured.")
endif()
  
message("libCGAL_ImageIO is configured")

if(NOT CGAL_HEADER_ONLY AND NOT ZLIB_FOUND)
  message( STATUS "NOTICE: libCGAL_ImageIO needs ZLib to read compressed files. That feature will not be activated.")
endif()

