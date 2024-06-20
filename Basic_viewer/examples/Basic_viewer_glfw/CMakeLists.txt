# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Basic_Viewer_GLFW)

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Eigen3 3.1.0)
include(CGAL_Eigen3_support)

set(GLFW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../include/CGAL/GLFW/vendor/glfw)
set(GLAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../include/CGAL/GLFW/vendor/glad)

include_directories(
    ${GLAD_DIR}/include
    ${GLFW_DIR}/include
    ${GLFW_DIR}/deps/ # for stb_image
)

add_subdirectory(${GLFW_DIR} ${CMAKE_BINARY_DIR}/glfw)

file(GLOB GLAD_SOURCES ${GLAD_DIR}/src/glad.c)
add_library(glad STATIC ${GLAD_SOURCES})

# file(GLOB GLFW_SOURCES ${GLFW_DIR}/src/*.c)
# add_library(glfw STATIC ${GLFW_SOURCES})

add_executable(${PROJECT_NAME} main.cpp)
add_executable(draw_mesh_and_points draw_mesh_and_points.cpp)
add_executable(draw_surface_mesh draw_surface_mesh.cpp)
add_executable(draw_surface_mesh_height draw_surface_mesh_height.cpp)
add_executable(screenshot screenshot.cpp)

target_link_libraries(${PROJECT_NAME} CGAL::CGAL glfw glad)
target_link_libraries(draw_mesh_and_points CGAL::CGAL glfw glad)
target_link_libraries(draw_surface_mesh CGAL::CGAL glfw glad)
target_link_libraries(draw_surface_mesh_height CGAL::CGAL glfw glad)
target_link_libraries(screenshot CGAL::CGAL glfw glad)

if(TARGET CGAL::Eigen3_support)
  target_link_libraries(${PROJECT_NAME} CGAL::Eigen3_support)
  target_link_libraries(draw_surface_mesh_height CGAL::Eigen3_support)
  target_link_libraries(draw_mesh_and_points CGAL::Eigen3_support)
endif()

add_definitions(
  -DCGAL_USE_BASIC_VIEWER_GLFW
  -DGLFW_INCLUDE_NONE
  -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"
)
