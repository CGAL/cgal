# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_sweep_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

set(CGAL_SEGMENT_TRAITS 1)
set(CGAL_POLYLINE_TRAITS 11)
set(CGAL_CONIC_TRAITS 21)
set(CGAL_POLYCONIC_TRAITS 22)

set(TRAP 1) # Trapezoidal decomposition
set(NAIVE 2)
set(WALK 3)

function(compile name source_file point_location traits)
  #    message(compile ${source_file})
  add_executable(${name} ${source_file})
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS ${name})
  if(TARGET CGAL::CGAL_Core)
    target_link_libraries(${name} CGAL::CGAL CGAL::CGAL_Core)
  else()
    target_link_libraries(${name} CGAL::CGAL)
  endif()
  target_compile_definitions(
    ${name} PRIVATE -DCGAL_ARR_TEST_POINT_LOCATION=${point_location}
                    -DCGAL_ARR_TEST_TRAITS=${traits})
endfunction()

function(compile_and_run_sweep name source_file point_location traits data_set)
  compile(${ARGV})
  if(NOT CGAL_ENABLE_TESTING)
    return()
  endif()
  cgal_add_compilation_test(${name})
  file(
    GLOB files
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/${data_set}/*.txt")
  foreach(file ${files})
    # message("test ${source_file} ${file}")
    string(MAKE_C_IDENTIFIER "${name}  ${file}" test_name)
    # message(" --> ${test_name}")
    cgal_add_test(${name} TEST_NAME ${test_name} ARGUMENTS ${file})
  endforeach()
endfunction()

compile_and_run_sweep(test_sweep test_sweep.cpp ${NAIVE}
                      ${CGAL_SEGMENT_TRAITS} "data/segments_tight")
compile_and_run_sweep(test_sweep_conic test_sweep.cpp ${NAIVE}
                      ${CGAL_CONIC_TRAITS} "data/conics")
compile_and_run_sweep(test_sweep_polyline test_sweep.cpp ${NAIVE}
                      ${CGAL_POLYLINE_TRAITS} "data/polylines")
