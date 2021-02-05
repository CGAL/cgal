# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Triangulation_2_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

set(CMAKE_AUTOMOC TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)
find_package(Qt5 QUIET COMPONENTS Widgets)

if(NOT CGAL_Qt5_FOUND
   OR NOT Qt5_FOUND)
  message(
    STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")
  return()
endif()

add_definitions(-DQT_NO_KEYWORDS)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#--------------------------------
# The "constrained Delaunay" demo: Constrained_Delaunay_triangulation_2
#--------------------------------

qt5_add_resources(CD_RES_FILE Constrained_Delaunay_triangulation_2.qrc)
qt5_wrap_ui(CD_UI_FILES Constrained_Delaunay_triangulation_2.ui)
# The executable itself.
add_executable(
  Constrained_Delaunay_triangulation_2
  Constrained_Delaunay_triangulation_2.cpp ${CD_UI_FILES}
  ${CGAL_Qt5_RESOURCE_FILES} ${CD_RES_FILE})
target_link_libraries(Constrained_Delaunay_triangulation_2
                      PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)
target_include_directories(Constrained_Delaunay_triangulation_2
                           PRIVATE ./include)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Constrained_Delaunay_triangulation_2)

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------
qt5_wrap_ui(D_UI_FILES Delaunay_triangulation_2.ui)
qt5_add_resources(D_RES_FILE Delaunay_triangulation_2.qrc)
# The executable itself.
add_executable(
  Delaunay_triangulation_2 Delaunay_triangulation_2.cpp ${D_UI_FILES}
                           ${CGAL_Qt5_RESOURCE_FILES} ${D_RES_FILE})
target_link_libraries(Delaunay_triangulation_2
                      PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Widgets)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_2)

#--------------------------------
# The "Regular" demo: Regular_triangulation_2
#--------------------------------

# The executable itself.
qt5_add_resources(R_RES_FILE Regular_triangulation_2.qrc)
qt5_wrap_ui(R_UI_FILES Regular_triangulation_2.ui)
add_executable(
  Regular_triangulation_2 Regular_triangulation_2.cpp ${R_UI_FILES}
                          ${CGAL_Qt5_RESOURCE_FILES} ${R_RES_FILE})
target_link_libraries(Regular_triangulation_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                                      Qt5::Widgets)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Regular_triangulation_2)

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
foreach(target Constrained_Delaunay_triangulation_2 Delaunay_triangulation_2
               Regular_triangulation_2)
  cgal_add_compilation_test(${target})
endforeach()
