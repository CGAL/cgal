message(STATUS "Using CGAL Lab from: ${CGAL_LAB_DEMO_DIR}")

list(INSERT CMAKE_MODULE_PATH 0 "${CGAL_LAB_DEMO_CMAKE_MODULE_PATH}")

message(STATUS "Now CMAKE_MODULE_PATH is: ${CMAKE_MODULE_PATH}")

include( CGALlab_macros )
include( cgal_lab_targets )


include_directories( ${CGAL_LAB_DEMO_HEADERS_DIRS} )
add_definitions(${CGAL_LAB_DEMO_DEFINITIONS})
