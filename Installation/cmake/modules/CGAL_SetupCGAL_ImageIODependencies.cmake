if(CGAL_SetupCGAL_ImageIODependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_ImageIODependencies_included TRUE)

find_package( ZLIB )

define_property(TARGET PROPERTY CGAL_TARGET_USES_ZLIB
  BRIEF_DOCS "Tells if the target uses ZLIB as a dependency"
  FULL_DOCS "Tells if the target uses ZLIB as a dependency")

if(ZLIB_FOUND)
  cache_set(CGAL_ImageIO_USE_ZLIB "ON")
endif(ZLIB_FOUND)

set( CGAL_ImageIO_BASENAME CGAL_ImageIO)

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

function(CGAL_setup_CGAL_ImageIO_dependencies target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()

  target_include_directories( CGAL_ImageIO SYSTEM ${keyword} ${ZLIB_INCLUDE_DIRS})
  target_link_libraries( CGAL_ImageIO ${keyword} CGAL::CGAL)
  target_link_libraries( CGAL_ImageIO ${keyword} ${ZLIB_LIBRARIES})

  if(ZLIB_FOUND)
    target_compile_definitions( CGAL_ImageIO ${keyword} "-DCGAL_USE_ZLIB")
    if(NOT ARGV1 STREQUAL INTERFACE)
      set_target_properties(CGAL_ImageIO PROPERTIES CGAL_TARGET_USES_ZLIB TRUE)
    endif()
  else()
    message( STATUS "NOTICE: libCGAL_ImageIO needs ZLib to read compressed files. That feature will not be activated.")
  endif()
endfunction()
