#.rst:
# CGAL_SetupCGAL_ImageIODependencies
# ----------------------------------
#
# The module searches for the dependencies of the `CGAL_ImageIO` library:
#   - the `Zlib` library (optional)
#
# by calling
#
# .. code-block:: cmake
#
#    find_package(ZLIB)
#
# and defines the variable :variable:`CGAL_ImageIO_FOUND` and the function
# :command:`CGAL_setup_CGAL_ImageIO_dependencies`.
#
if(CGAL_SetupCGAL_ImageIODependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_ImageIODependencies_included TRUE)

#.rst:
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`FindZLIB`
find_package( ZLIB )

define_property(TARGET PROPERTY CGAL_TARGET_USES_ZLIB
  BRIEF_DOCS "Tells if the target uses ZLIB as a dependency"
  FULL_DOCS "Tells if the target uses ZLIB as a dependency")

#.rst:
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_ImageIO_USE_ZLIB
#
#    Set to `TRUE` if `CGAL_ImageIO` was compiled with `Zlib`.
#
# .. variable:: CGAL_ImageIO_FOUND
#
#    Always set to `TRUE`.

if(ZLIB_FOUND)
  set(CGAL_ImageIO_USE_ZLIB ON CACHE BOOL "CGAL_ImageIO uses ZLIB")
endif(ZLIB_FOUND)

set( CGAL_ImageIO_BASENAME CGAL_ImageIO)

set(CGAL_ImageIO_FOUND TRUE)
set_property(GLOBAL PROPERTY CGAL_ImageIO_FOUND TRUE)

function(CGAL_setup_CGAL_ImageIO_dependencies target)

  target_link_libraries( CGAL_ImageIO INTERFACE CGAL::CGAL)

  if(ZLIB_FOUND)
    target_include_directories( CGAL_ImageIO SYSTEM INTERFACE ${ZLIB_INCLUDE_DIRS})
    target_link_libraries( CGAL_ImageIO INTERFACE ${ZLIB_LIBRARIES})

    target_compile_definitions( CGAL_ImageIO INTERFACE CGAL_USE_ZLIB=1)
  else()
    message( STATUS "NOTICE: libCGAL_ImageIO needs ZLib to read compressed files. That feature will not be activated.")
  endif()
endfunction()
