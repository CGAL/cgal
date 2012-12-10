project(LCC_performance)

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

include_directories(BEFORE "../../include")

add_executable(lcc_performance performance.h lcc_performance.h lcc_performance.cpp)
target_link_libraries(lcc_performance ${CGAL_LIBRARIES})
