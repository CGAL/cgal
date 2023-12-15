# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Snap_rounding_2_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program(test_snap_rounding_2.cpp NO_TESTING)

function(add_Snap_rounding_tests name)
  set(data_dir "data")
  file(
    GLOB files
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${data_dir}
    "${CMAKE_CURRENT_SOURCE_DIR}/${data_dir}/data*")
  foreach(file ${files})
    string(MAKE_C_IDENTIFIER "${name}  ${data_dir}/${file}" test_name)
    cgal_add_test(${name} TEST_NAME ${test_name} ARGUMENTS ${data_dir}/${file}
                  ${data_dir}/output_${file})
    add_test(NAME "${test_name}_compare_results"
             COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
                     ${data_dir}/output_${file} ${data_dir}/gold_${file})
    cgal_setup_test_properties("${test_name}_compare_results")
    set_property(
      TEST "${test_name}_compare_results"
      APPEND
      PROPERTY DEPENDS "${test_name}")
    if(POLICY CMP0066) # CMake 3.7 or later
      set_property(TEST "${test_name}" PROPERTY FIXTURES_SETUP
                                                "${test_name}_compare_results")
      set_property(TEST "${test_name}_compare_results" APPEND
                   PROPERTY FIXTURES_REQUIRED "${test_name}_compare_results")
    endif()
  endforeach()
endfunction(add_Snap_rounding_tests)

if(CGAL_ENABLE_TESTING)
  cgal_add_compilation_test(test_snap_rounding_2)
  add_snap_rounding_tests(test_snap_rounding_2)
endif()
