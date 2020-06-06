# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( AABB_tree )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )
FIND_PACKAGE(embree 3.0 REQUIRED)

include_directories(include)

ADD_EXECUTABLE(AabbTest AabbTest.cpp include/RaysGenerate.h)
TARGET_LINK_LIBRARIES(AabbTest CGAL)

ADD_EXECUTABLE(EmbreeTest EmbreeTest.cpp include/RaysGenerate.h)
TARGET_LINK_LIBRARIES(EmbreeTest embree)


