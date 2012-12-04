option( WITH_GMP "Use the GMP number types if available." ON )
if ( WITH_GMP )
  include(CGAL_SetupGMPXX)
  include(CGAL_SetupGMP)
endif( WITH_GMP )

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
  message( STATUS "CGAL_Core needs GMP, cannot be configured.")
endif( NOT GMP_FOUND )

option ( WITH_LEDA "Use the LEDA number types if available." OFF )
if ( WITH_LEDA )
  include(CGAL_SetupLEDA)
endif( WITH_LEDA )

option ( WITH_MPFI "Use MPFI if available." OFF )
if ( WITH_MPFI )
  include(CGAL_SetupMPFI)
endif( WITH_MPFI )

option ( WITH_RS "Use RS if available." OFF )
if ( WITH_RS )
  include(CGAL_SetupRS)
endif( WITH_RS )

option ( WITH_NTL "Use NTL if available." OFF )
if ( WITH_NTL )
  include(CGAL_SetupNTL)
endif( WITH_NTL )

include(CGAL_SetupBoost)
