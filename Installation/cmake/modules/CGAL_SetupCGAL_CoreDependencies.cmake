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
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`CGAL_SetupGMP`
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_Core_FOUND
#
#    Set to `TRUE` if the dependencies of `CGAL_Core` were found.

if(NOT CGAL_DISABLE_GMP)
  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupGMP.cmake)
  if(GMP_FOUND)
    set(CGAL_Core_FOUND TRUE)
    set_property(GLOBAL PROPERTY CGAL_Core_FOUND TRUE)
  endif()
endif()

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
  use_CGAL_GMP_support(CGAL_Core INTERFACE)
  target_compile_definitions(${target} INTERFACE CGAL_USE_CORE=1)
  target_link_libraries( CGAL_Core INTERFACE CGAL::CGAL )

endfunction()
