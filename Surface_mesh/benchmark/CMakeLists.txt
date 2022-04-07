cmake_minimum_required(VERSION 3.1...3.20)
project(Surface_mesh_performance)

find_package(CGAL REQUIRED)

include_directories(BEFORE "../include")

# For profilling with gprof
#add_definitions("-pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# add_definitions("-g")
add_definitions("-std=c++1y")

# Polyhedron
add_executable(polyhedron_performance performance_2.h polyhedron_performance.h
                                      polyhedron_performance.cpp)
target_link_libraries(polyhedron_performance ${CGAL_LIBRARIES})

# LCC_2
add_executable(lcc_performance_2 performance_2.h lcc_performance_2.h
                                 lcc_performance_2.cpp)
target_link_libraries(lcc_performance_2 ${CGAL_LIBRARIES})

# Surface_mesh
add_executable(
  surface_mesh_performance performance_2.h surface_mesh_performance.h
                           surface_mesh_performance.cpp)

# Performance_2
add_executable(
  performance_2 performance_2.cpp performance_2.h polyhedron_performance.h
                surface_mesh_performance.h lcc_performance_2.h)

target_link_libraries(performance_2 ${CGAL_LIBRARIES})

create_single_source_cgal_program("sm_sms.cpp")
create_single_source_cgal_program("poly_sms.cpp")
