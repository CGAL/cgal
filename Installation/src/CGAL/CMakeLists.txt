message("Configuring libCGAL")

use_essential_libs()

include_directories (SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

link_directories(${CGAL_3RD_PARTY_LIBRARIES_DIRS})

collect_cgal_library(CGAL "")

if(CGAL_HEADER_ONLY)
  target_link_libraries(CGAL INTERFACE ${CGAL_3RD_PARTY_LIBRARIES})
else()
  target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})
endif()

message("libCGAL is configured")
