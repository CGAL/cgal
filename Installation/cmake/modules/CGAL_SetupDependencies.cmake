include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Macros.cmake)

if( (GMP_FOUND AND NOT MPFR_FOUND) OR (NOT GMP_FOUND AND MPFR_FOUND) )
  message( FATAL_ERROR "CGAL needs for its full functionality both GMP and MPFR.")
endif()

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
  message( STATUS "CGAL_Core needs GMP, cannot be configured.")
endif( NOT GMP_FOUND )

# finally setup Boost
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupBoost.cmake)
