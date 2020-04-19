# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Three_Demo )


set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
    include(${CGAL_USE_FILE})
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

