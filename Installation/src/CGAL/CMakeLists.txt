include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

link_directories(${CGAL_3RD_PARTY_LIBRARIES_DIRS})

build_cgal_library(CGAL CGAL "") 

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

