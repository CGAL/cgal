#.rst:
# CGAL_SetupGMP
# -------------
#
# The module searchs for the `GMP` and `MPFR` headers and libraries,
# by calling
#
# .. code-block:: cmake
#
#    find_package(GMP)
#    find_package(MPFR)
#
# and defines the function :command:`use_CGAL_GMP_support`.

if(CGAL_SetupGMP_included OR CGAL_DISABLE_GMP)
  return()
endif()
set(CGAL_SetupGMP_included TRUE)

find_package(GMP)
find_package(MPFR)

#.rst:
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: use_CGAL_GMP_support
#
#    Link the target with the `GMP` and `MPFR` libraries::
#
#      use_CGAL_GMP_support( target [INTERFACE] )
#
#    If the option ``INTERFACE`` is passed, the dependencies are
#    added using :command:`target_link_libraries` with the ``INTERFACE``
#    keyword, or ``PUBLIC`` otherwise.

function(use_CGAL_GMP_support target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(NOT GMP_FOUND OR NOT MPFR_FOUND)
    message(FATAL_ERROR "CGAL requires GMP and MPFR.")
    return()
  endif()

  target_include_directories(${target} SYSTEM ${keyword} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
  target_link_libraries(${target} ${keyword} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
endfunction()
