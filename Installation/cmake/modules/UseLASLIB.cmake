# This module setups the compiler for using LASlib library.
# It assumes that find_package(LASLIB) was already called.

add_definitions(-DCGAL_LINKED_WITH_LASLIB)
