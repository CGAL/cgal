message(STATUS "Configure libCGAL")

use_lib( GMP "###")
use_lib( MPFR "###")
if (FOUND_GMPXX) 
  use_lib( GMPXX "###")
endif()

include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

link_directories(${CGAL_3RD_PARTY_LIBRARIES_DIRS})

collect_cgal_library(CGAL CGAL "") 

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

message(STATUS "libCGAL is configured")
