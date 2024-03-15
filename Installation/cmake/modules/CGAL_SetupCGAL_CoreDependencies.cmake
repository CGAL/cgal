#.rst:
# CGAL_SetupCGAL_CoreDependencies
# -------------------------------
#
# The module searches for the dependencies of the `CGAL_Core` library:
#   - the `GMP/MPFR` couple,
#
# and defines the variable :variable:`CGAL_Core_FOUND` and the function
# :command:`CGAL_setup_CGAL_Core_dependencies`.
#
# Module Input Variables
# ^^^^^^^^^^^^^^^^^^^^^^
# - :variable:`CGAL_DISABLE_GMP`

if(CGAL_SetupCGAL_CoreDependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_CoreDependencies_included TRUE)

#.rst:
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_Core_FOUND
#
#    Set to `TRUE` if the dependencies of `CGAL_Core` were found.


# always found as it requires the minimal version of boost required by CGAL
set(CGAL_Core_FOUND TRUE)
set_property(GLOBAL PROPERTY CGAL_Core_FOUND TRUE)

#.rst:
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: CGAL_setup_CGAL_Core_dependencies
#
#   Link the target with the dependencies of `CGAL_Core`::
#
#     CGAL_setup_CGAL_Core_dependencies( target)
#
#   The dependencies are
#   added using :command:`target_link_libraries` with the ``INTERFACE``
#   keyword.
#

function(CGAL_setup_CGAL_Core_dependencies target)
  find_package( Boost 1.72 REQUIRED )
  if (!CGAL_DISABLE_GMP AND GMP_FOUND)
    use_CGAL_GMP_support(CGAL_Core INTERFACE)
  endif()
  target_compile_definitions(${target} INTERFACE CGAL_USE_CORE=1)
  target_link_libraries( CGAL_Core INTERFACE CGAL::CGAL )
endfunction()
