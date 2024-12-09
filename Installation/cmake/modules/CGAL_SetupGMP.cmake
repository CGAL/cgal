#.rst:
# CGAL_SetupGMP
# -------------
#
# The module searches for the `GMP` and `MPFR` headers and libraries,
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

# Locally setting of CMAKE_MODULE_PATH, not exported to parent scope.
# That is required to find the FindGMP and FindMPFR modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_MODULES_DIR})

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(GMPXX QUIET)

if(NOT GMPXX_FOUND)
  option(CGAL_WITH_GMPXX "Use CGAL with GMPXX: use C++ classes of GNU MP instead of CGAL wrappers" OFF)
else()
  option(CGAL_WITH_GMPXX "Use CGAL with GMPXX: use C++ classes of GNU MP instead of CGAL wrappers" ON)
endif()

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
  if(NOT GMP_FOUND OR NOT MPFR_FOUND)
    message(FATAL_ERROR "CGAL requires GMP and MPFR.")
    return()
  endif()

  if(NOT GMP_INCLUDE_DIR STREQUAL "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include")
    target_include_directories(${target} SYSTEM INTERFACE ${GMP_INCLUDE_DIR})
  else()
    target_include_directories(${target} SYSTEM INTERFACE
      $<BUILD_INTERFACE:${GMP_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:include>)
  endif()
  if(NOT MPFR_INCLUDE_DIR STREQUAL "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include")
    target_include_directories(${target} SYSTEM INTERFACE ${MPFR_INCLUDE_DIR})
  else()
    target_include_directories(${target} SYSTEM INTERFACE
      $<BUILD_INTERFACE:${MPFR_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:include>)
  endif()
  if(WITH_GMPXX OR CGAL_WITH_GMPXX)
    target_include_directories(${target} SYSTEM INTERFACE  ${GMPXX_INCLUDE_DIR})
    target_link_libraries(${target}  INTERFACE ${GMPXX_LIBRARIES})
    target_compile_definitions(${target} INTERFACE CGAL_USE_GMPXX=1)
  endif()
  target_link_libraries(${target} INTERFACE ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
endfunction()
