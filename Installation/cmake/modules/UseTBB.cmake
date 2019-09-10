# This module setups the compiler for using TBB library.
# It assumes that find_package(TBB) was already called.

include_directories ( ${TBB_INCLUDE_DIRS} )
link_directories( ${TBB_LIBRARY_DIRS} )
add_definitions( -DNOMINMAX -DCGAL_LINKED_WITH_TBB )

message(DEPRECATION "This file UseTBB.cmake is deprecated, and the function `CGAL_target_use_TBB` from CGAL_target_use_TBB.cmake should be used instead.")
