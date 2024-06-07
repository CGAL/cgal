cmake_minimum_required(VERSION 3.1...3.23)
project(Mesh_3)


if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/${PROJECT_NAME}/CMakeLists.txt")
  option(CGAL_${PROJECT_NAME}_BENCHMARK "Add benchmark to project" OFF)
  if(CGAL_${PROJECT_NAME}_BENCHMARK)
    add_subdirectory("benchmark/${PROJECT_NAME}")
  endif()
endif()

add_subdirectory("examples/${PROJECT_NAME}")
add_subdirectory("test/${PROJECT_NAME}")
