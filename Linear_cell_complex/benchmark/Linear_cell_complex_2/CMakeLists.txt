project(LCC_performance_2)

cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

add_subdirectory(surface_mesh)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

include_directories(BEFORE "../../include")
include_directories(BEFORE "./surface_mesh")
include_directories(BEFORE "./cgogn" "./cgogn/include")
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cgogn/lib/Release)
include_directories(BEFORE "/usr/include/libxml2/")

# For profilling with gprof
#add_definitions("-pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# add_definitions("-g")

# OpenMesh
set(OPENMESH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmesh")
find_package(OpenMesh REQUIRED)
include_directories(${OPENMESH_INCLUDE_DIR})
link_directories(${OPENMESH_LIBRARY_DIR})

# Polyhedron
add_executable(polyhedron_performance performance_2.h polyhedron_performance.h polyhedron_performance.cpp)
target_link_libraries(polyhedron_performance ${CGAL_LIBRARIES})

# LCC_2
add_executable(lcc_performance_2 performance_2.h lcc_performance_2.h lcc_performance_2.cpp)
target_link_libraries(lcc_performance_2 ${CGAL_LIBRARIES})

# Surface_mesh
add_executable(surface_mesh_performance performance_2.h surface_mesh_performance.h surface_mesh_performance.cpp)
target_link_libraries(surface_mesh_performance surface_mesh)

# Open_mesh
add_executable(openmesh_performance performance_2.h openmesh_performance.h openmesh_performance.cpp)
target_link_libraries(openmesh_performance ${OPENMESH_LIBRARIES})

# CGoGN
find_package(Qt REQUIRED)
SET(QT_USE_QTSVG TRUE)
SET(QT_USE_QTXML TRUE )
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
SET (CGoGN_EXT_INCLUDES ${CGoGN_EXT_INCLUDES} ${QT_INCLUDE_DIR} ${QGLVIEWER_INCLUDE_DIR} )
SET (CGoGN_EXT_LIBS ${CGoGN_EXT_LIBS} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )

add_executable(cgogn_performance_2 performance_2.h cgogn_performance_2.h cgogn_performance_2.cpp)
target_link_libraries(cgogn_performance_2 algo assimp container nl topology utils Zinri z ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )

# Performance_2
add_executable(performance_2
  performance_2.cpp
  performance_2.h
  polyhedron_performance.h
  openmesh_performance.h
  surface_mesh_performance.h
  cgogn_performance_2.h
  lcc_performance_2.h
)

target_link_libraries(performance_2
  ${CGAL_LIBRARIES}
  surface_mesh
  algo assimp container nl topology utils Zinri z ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES}
  ${OPENMESH_LIBRARIES}
)
