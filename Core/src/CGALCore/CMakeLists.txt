message(STATUS "Configure libCGALCore")
 
use_lib( GMP "###")
use_lib( MPFR "###")
if (FOUND_GMPXX) 
  use_lib( GMPXX "###")
endif()

include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

collect_cgal_library(CGAL_Core CGALCore "")
 
target_link_libraries( CGAL_Core ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} )

add_dependencies( CGAL_Core CGAL )

message(STATUS "libCGALCore is configured")

