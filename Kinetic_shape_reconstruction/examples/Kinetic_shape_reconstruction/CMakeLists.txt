# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(KSR_Examples)

cmake_minimum_required(VERSION 3.1...3.15)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)
  message(STATUS "Found CGAL")

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    message(STATUS "Found Boost")

    set(targets
    # kinetic_2d_example
    kinetic_precomputed_shapes_example
    # kinetic_random_shapes_example
    )

    set(project_linked_libraries)
    set(project_compilation_definitions)

    foreach(target ${targets})
      create_single_source_cgal_program("${target}.cpp")
      if(TARGET ${target})
        target_link_libraries(${target} PUBLIC ${project_linked_libraries})
        target_compile_definitions(${target} PUBLIC ${project_compilation_definitions})
      endif()
    endforeach()
  else()
    message(ERROR "This program requires the Boost library, and will not be compiled.")
  endif()
else()
  message(ERROR "This program requires the CGAL library, and will not be compiled.")
endif()
