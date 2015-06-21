message(STATUS "Using CGAL_polyhedron_demo from: ${CGAL_POLYHEDRON_DEMO_DIR}")

list(INSERT CMAKE_MODULE_PATH 0 "${CGAL_POLYHEDRON_DEMO_CMAKE_MODULE_PATH}")

message(STATUS "Now CMAKE_MODULE_PATH is: ${CMAKE_MODULE_PATH}")

include( polyhedron_demo_macros )
include( polyhedron_demo_targets )


include_directories( ${CGAL_POLYHEDRON_DEMO_HEADERS_DIRS} )
add_definitions(${CGAL_POLYHEDRON_DEMO_DEFINITIONS})
