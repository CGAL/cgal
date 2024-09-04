cmake_minimum_required(VERSION 3.12...3.29)
project(LCC_performance_2)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

add_subdirectory(surface_mesh)

find_package(CGAL REQUIRED)

include_directories(BEFORE "./surface_mesh")
include_directories(BEFORE "./cgogn" "./cgogn/include")
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/cgogn/lib/Release)
include_directories(BEFORE "/usr/include/libxml2/")

# For profilling with gprof
#add_definitions("-pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# add_definitions("-g")

# OpenMesh
find_package(OpenMesh REQUIRED)
include(CGAL_OpenMesh_support)

# Polyhedron
add_executable(polyhedron_performance performance_2.h polyhedron_performance.h
                                      polyhedron_performance.cpp)
target_link_libraries(polyhedron_performance CGAL::CGAL)

# LCC_2
add_executable(lcc_performance_2 performance_2.h lcc_performance_2.h
                                 lcc_performance_2.cpp)
target_link_libraries(lcc_performance_2 CGAL::CGAL)

# Surface_mesh
add_executable(
  surface_mesh_performance performance_2.h surface_mesh_performance.h
                           surface_mesh_performance.cpp)
target_link_libraries(surface_mesh_performance surface_mesh)

# Open_mesh
add_executable(openmesh_performance performance_2.h openmesh_performance.h
                                    openmesh_performance.cpp)
target_link_libraries(openmesh_performance CGAL::OpenMesh_support)

# CGoGN
find_package(Qt REQUIRED)
set(QT_USE_QTSVG TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
set(CGoGN_EXT_INCLUDES ${CGoGN_EXT_INCLUDES} ${QT_INCLUDE_DIR})
set(CGoGN_EXT_LIBS ${CGoGN_EXT_LIBS} ${QT_LIBRARIES})

add_executable(cgogn_performance_2 performance_2.h cgogn_performance_2.h
                                   cgogn_performance_2.cpp)
target_link_libraries(
  cgogn_performance_2
  algo
  assimp
  container
  nl
  topology
  utils
  Zinri
  z
  ${QT_LIBRARIES})

# Performance_2
add_executable(
  performance_2
  performance_2.cpp
  performance_2.h
  polyhedron_performance.h
  openmesh_performance.h
  surface_mesh_performance.h
  cgogn_performance_2.h
  lcc_performance_2.h)

target_link_libraries(
  performance_2
  CGAL::CGAL
  surface_mesh
  algo
  assimp
  container
  nl
  topology
  utils
  Zinri
  z
  ${QT_LIBRARIES}
  CGAL::OpenMesh_support)
