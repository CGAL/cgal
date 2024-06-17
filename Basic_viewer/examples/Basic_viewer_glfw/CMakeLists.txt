# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(GLFW_Basicv)

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Eigen3 3.1.0)
include(CGAL_Eigen3_support)

include_directories (
    .
    ../../include/CGAL/GLFW/vendor/glfw/deps/
    ../../include/CGAL/GLFW/vendor/glad/include/
    ../../include/CGAL/GLFW/vendor/glfw/include/
)


file (GLOB         VENDORS_SOURCES ../../include/CGAL/GLFW/vendor/glad/src/glad.c)

file (GLOB_RECURSE PROJECT_HEADERS *.hpp
                                   *.h)
file (GLOB_RECURSE PROJECT_SOURCES *.cpp
                                   *.cxx
                                   *.cc
                                   *.c)

add_executable (${PROJECT_NAME} ${VENDORS_SOURCES} "main.cpp")
add_executable ("draw_mesh_and_points" ${VENDORS_SOURCES} "draw_mesh_and_points.cpp")
add_executable ("draw_surface_mesh" ${VENDORS_SOURCES} "draw_surface_mesh.cpp")
add_executable ("draw_surface_mesh_height" ${VENDORS_SOURCES} "draw_surface_mesh_height.cpp")
add_executable ("screenshot" ${VENDORS_SOURCES} "screenshot.cpp")

target_link_libraries(${PROJECT_NAME} glfw CGAL::CGAL)
target_link_libraries(draw_mesh_and_points glfw CGAL::CGAL)
target_link_libraries(draw_surface_mesh glfw CGAL::CGAL)
target_link_libraries(draw_surface_mesh_height glfw CGAL::CGAL)
target_link_libraries(screenshot glfw CGAL::CGAL)

if(TARGET CGAL::Eigen3_support)
  target_link_libraries(GLFW_Basicv CGAL::Eigen3_support)
  target_link_libraries(draw_surface_mesh_height CGAL::Eigen3_support)
  target_link_libraries(draw_mesh_and_points CGAL::Eigen3_support)
endif()

add_definitions (
  -DCGAL_USE_BASIC_VIEWER_GLFW
  -DGLFW_INCLUDE_NONE
  -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"
)
