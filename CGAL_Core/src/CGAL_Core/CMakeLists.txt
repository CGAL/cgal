message("Configuring libCGAL_Core")
 
use_essential_libs()

include_directories (SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )

collect_cgal_library(CGAL_Core "")
 
target_link_libraries( CGAL_Core ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} )

add_dependencies( CGAL_Core CGAL )

message("libCGAL_Core is configured")

