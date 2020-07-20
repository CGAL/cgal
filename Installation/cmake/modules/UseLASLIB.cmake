# This module setups the compiler for using LASlib library.
# It assumes that find_package(LASLIB) was already called.

add_definitions(-DCGAL_LINKED_WITH_LASLIB)

message(DEPRECATION "This file UseLASLIB.cmake is deprecated, and the imported target `CGAL::TBB_support` from CGAL_LASLIB_support.cmake should be used instead.")
