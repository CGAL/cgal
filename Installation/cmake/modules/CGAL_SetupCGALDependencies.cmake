#.rst:
# CGAL_SetupCGALDependencies
# --------------------------
#
# The module searchs for the dependencies of the CGAL library:
#   - the `GMP/MPFR` couple,
#   - `LEDA` (optional)
#   - the `Boost` libraries (mostly the header-only libraries)
#
# and defines the variable :variable:`CGAL_FOUND` and the function
# :command:`CGAL_setup_CGAL_dependencies`.
#
# Module Input Variables
# ^^^^^^^^^^^^^^^^^^^^^^
# .. variable:: CGAL_DISABLE_GMP
#
#    If set, the `GMP` library will not be used. If
#    :variable:`WITH_LEDA` is not used either, a efficient exact
#    number types are used by CGAL kernels for exact computation.
#
# .. variable:: WITH_LEDA
#
#    If set, the `LEDA` library will be searched and used to provide
#    the exact number types used by CGAL kernels.
#
# .. variable:: CGAL_HEADER_ONLY
#
#    Set this variable if you are using the CGAL libraries as
#    header-only libraries.
#
if(CGAL_SetupCGALDependencies_included)
  return()
endif()
set(CGAL_SetupCGALDependencies_included TRUE)

#.rst:
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`CGAL_SetupGMP`
if(NOT CGAL_DISABLE_GMP)
  include(CGAL_SetupGMP)
endif()

#.rst:
#   - :module:`CGAL_SetupLEDA`
if(WITH_LEDA)
  include(CGAL_SetupLEDA)
endif()

#.rst:
#   - :module:`CGAL_SetupBoost`
include(CGAL_SetupBoost)

#.rst:
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_FOUND
#
#    Set to `TRUE` if the dependencies of CGAL were found.
if(Boost_FOUND)
  set(CGAL_FOUND TRUE)
endif()

#.rst:
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: CGAL_setup_CGAL_dependencies
#
#   Link the target with the dependencies of CGAL::
#
#     CGAL_setup_CGAL_dependencies( target [INTERFACE] )
#
#   If the option ``INTERFACE`` is passed, the dependencies are
#   added using :command:`target_link_libraries` with the ``INTERFACE``
#   keyword, or ``PUBLIC`` otherwise.
#
function(CGAL_setup_CGAL_dependencies target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(NOT CGAL_DISABLE_GMP)
    use_CGAL_GMP_support(${target} ${keyword})
    set(CGAL_USE_GMP  TRUE CACHE INTERNAL "CGAL library is configured to use GMP")
    set(CGAL_USE_MPFR TRUE CACHE INTERNAL "CGAL library is configured to use MPFR")
  endif()

  if(WITH_LEDA)
    use_CGAL_LEDA_support(${target} ${keyword})
  endif()
  
  if (CGAL_HEADER_ONLY)
    target_compile_definitions(${target} ${keyword} CGAL_HEADER_ONLY=1)
  endif()

  use_CGAL_Boost_support(${target} ${keyword})
  foreach(dir ${CGAL_INCLUDE_DIRS})
    target_include_directories(${target} ${keyword}
      $<BUILD_INTERFACE:${dir}>)
  endforeach()
  target_include_directories(${target} ${keyword}
    $<INSTALL_INTERFACE:include/CGAL>)
endfunction()
