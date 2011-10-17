include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

link_directories(${CGAL_3RD_PARTY_LIBRARIES_DIRS})

build_cgal_library(CGAL CGAL "") 

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

set( CGAL_BASENAME CGAL )

set_target_properties(CGAL PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

