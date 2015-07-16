project(LCC_performance_3)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

find_package(Boost 1.43.0)
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else()
   set (USE_IN_SOURCE_TREE_BOOST true)
endif()

add_subdirectory(openvolumemesh)

include_directories(BEFORE openvolumemesh/src)
include_directories(BEFORE "../../include")
include_directories(BEFORE "./cgogn" "./cgogn/include")
include_directories(BEFORE "/usr/include/libxml2/")
include_directories(BEFORE "/usr/include/eigen3/")

LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cgogn/lib/Release)

find_package(Qt REQUIRED)
SET(QT_USE_QTSVG TRUE)
SET(QT_USE_QTXML TRUE )
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

add_definitions(-DINCLUDE_TEMPLATES)

# Performance_3
add_executable(performance_3 performance_3.cpp)
add_dependencies(performance_3 OpenVolumeMesh)
target_link_libraries(performance_3 OpenVolumeMesh boost_timer boost_system ${CGAL_LIBRARIES}  algo assimp container nl topology utils Zinri z xml2 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${MAP_VIEWER_LIBRARIES})

# CGoGN
add_executable(cgogn_performance_3 performance_3.h cgogn_performance_3.h cgogn_performance_3.cpp)
target_link_libraries(cgogn_performance_3  algo assimp container nl topology utils Zinri z xml2 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )

# LCC_3
add_executable(lcc_performance_3 performance_3.h lcc_performance_3.h lcc_performance_3.cpp)
target_link_libraries(lcc_performance_3 ${CGAL_LIBRARIES} ${MAP_VIEWER_LIBRARIES})

# OpenVolumeMesh
add_executable(openvolumemesh_performance performance_3.h openvolumemesh_performance.h openvolumemesh_performance.cpp)
target_link_libraries(openvolumemesh_performance OpenVolumeMesh boost_timer boost_system)
