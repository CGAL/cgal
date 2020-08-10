# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project(Octree_Tests)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
if (CGAL_FOUND)

    # TODO: Add single-source programs

else ()
    message(WARNING
            "This program requires the CGAL library, and will not be compiled.")
endif ()
