# This CMake script is a translation of the bash script `./cgal_test_base`
# to the CMake language. It defines the targets to compile, as well as the
# tests to run with CTest.


# SET PARAMETERS FOR cgal_test

include(CGAL_add_test)

function(cgal_debug_message)
#  message(${ARGN})
endfunction()

set(ERRORFILE error.txt)
set(DO_RUN y)

set(FULL_ERROR_DESCRIPTION_FILE ProgramOutput.error.txt)

#---------------------------------------------------------------------#
#                    compile_and_run <target>
#---------------------------------------------------------------------#

# note that these values shloud match to the values in test_configuration.h file

set(CARTESIAN_KERNEL 0)
set(SIMPLE_CARTESIAN_KERNEL 1)
set(UNIVARIATE_ALGEBRAIC_KERNEL 2)

set(SEGMENT_GEOM_TRAITS 0)
set(NON_CACHING_SEGMENT_GEOM_TRAITS 1)
set(POLYLINE_GEOM_TRAITS 2)
set(NON_CACHING_POLYLINE_GEOM_TRAITS 3)
set(LINEAR_GEOM_TRAITS 4)
set(CORE_CONIC_GEOM_TRAITS 5)
set(LINE_ARC_GEOM_TRAITS 6)
set(CIRCULAR_ARC_GEOM_TRAITS 7)
set(CIRCULAR_LINE_ARC_GEOM_TRAITS 8)
set(CIRCLE_SEGMENT_GEOM_TRAITS 9)
set(BEZIER_GEOM_TRAITS 10)
set(GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS 11)
set(RATIONAL_ARC_GEOM_TRAITS 12)
set(ALGEBRAIC_GEOM_TRAITS 13)
set(POLYCURVE_CONIC_GEOM_TRAITS 14)
set(POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS 15)
set(POLYCURVE_BEZIER_GEOM_TRAITS 16)
set(FLAT_TORUS_GEOM_TRAITS 17)

set(PLANAR_BOUNDED_TOPOL_TRAITS 0)
set(PLANAR_UNBOUNDED_TOPOL_TRAITS 1)
set(SPHERICAL_TOPOL_TRAITS 2)

set(DOUBLE_NT 0)
set(MP_FLOAT_NT 1)
set(GMPZ_NT 2)
set(LEDA_RAT_NT 3)
set(QUOTIENT_MP_FLOAT_NT 4)
set(QUOTIENT_CGAL_GMPZ_NT 5)
set(CGAL_GMPQ_NT 6)
set(LAZY_LEDA_RAT_NT 7)
set(LAZY_CGAL_GMPQ_NT 8)
set(LAZY_QUOTIENT_MP_FLOAT_NT 9)
set(LEDA_REAL_NT 10)
set(CORE_EXPR_NT 11)
set(LAZY_GMPZ_NT 12)
set(LEDA_INT_NT 13)
set(CGAL_GMPZ_NT 14)
set(CORE_INT_NT 15)
set(CORE_RAT_NT 16)

if(CGAL_DISABLE_GMP)
  set(CGAL_DISABLE_GMP ON)
endif()

if(CGAL_DISABLE_GMP)
  message(STATUS "GMP is disable. Try to use LEDA instead.")
  set(GMPZ_NT ${LEDA_INT_NT})
  set(QUOTIENT_CGAL_GMPZ_NT ${LEDA_RAT_NT})
  set(CGAL_GMPQ_NT ${LEDA_RAT_NT})
  set(LAZY_CGAL_GMPQ_NT ${LAZY_LEDA_RAT_NT})
  set(LAZY_GMPZ_NT ${LAZY_LEDA_RAT_NT})
  set(CGAL_GMPZ_NT ${LEDA_INT_NT})
endif()

set(COMPARE 1)
set(VERTEX 2)
set(IS_VERTICAL 3)
set(COMPARE_Y_AT_X 4)
set(COMPARE_Y_AT_X_LEFT 5)
set(COMPARE_Y_AT_X_RIGHT 6)
set(MAKE_X_MONOTONE 7)
set(INTERSECT 8)
set(SPLIT 9)
set(ARE_MERGEABLE 10)
set(MERGE 11)
set(ASSERTIONS 12)
set(CONSTRUCTOR 13)
set(COMPARE_X_ON_BOUNDARY 16)
set(COMPARE_X_NEAR_BOUNDARY 17)
set(COMPARE_Y_NEAR_BOUNDARY 18)
set(PARAMETER_SPACE_X 19)
set(PARAMETER_SPACE_Y 20)
set(X_ON_IDENTIFICATION 21)
set(Y_ON_IDENTIFICATION 22)
set(IS_BOUNDED 23)
set(IS_IN_X_RANGE 24)
set(COMPARE_Y_POSITION 25)
set(IS_BETWEEN_CW 26)
set(COMPARE_CW_AROUND_POINT 27)
set(PUSH_BACK 28)
set(PUSH_FRONT 29)
set(NUMBER_OF_POINTS 32)
set(COMPARE_ENDPOINTS_XY 33)
set(CONSTRUCT_OPPOSITE 34)
set(TRIM 35)

#---------------------------------------------------------------------#
#                    configure
#---------------------------------------------------------------------#

function(configure type)
  cgal_debug_message(STATUS "---> cmake with -DCGAL_CXX_FLAGS:STRING=\"${TESTSUITE_CXXFLAGS}\"")
#  string(REPLACE "-DTEST_" "" suffix "${TESTSUITE_CXXFLAGS}")
#  string(REPLACE "KERNEL" "K" suffix "${suffix}")
#  string(REPLACE "GEOM_TRAITS" "GT" suffix "${suffix}")
#  string(REPLACE "-DCGAL_ARR_POINT_LOCATION" "POINT_LOCATION" suffix "${suffix}")
#  string(MAKE_C_IDENTIFIER "${suffix}" suffix)
#  cgal_debug_message(STATUS "       suffix: ${suffix}")
#  set(suffix ${suffix} PARENT_SCOPE)
  set(suffix ${type} PARENT_SCOPE)
endfunction()

function(compile_test_with_flags name type flags)
  cgal_debug_message(STATUS "# compile_test_with_flags(${name} ${type} \"${flags}\" ${ARGN})")
  set(TESTSUITE_CXXFLAGS "${flags}")
  set(TESTSUITE_CXXFLAGS "${flags}" PARENT_SCOPE)
#  message("   successful configuration")
#  message("   successful compilation of ${name} ${type}")
  if(ARGV3)
    configure(${ARGV3})
  else()
    configure(${type})
  endif()
  set(suffix ${suffix} PARENT_SCOPE)
#  set(suffix ${type})
  cgal_arr_2_add_target(${name} ${name}.cpp)
endfunction()

function(cgal_arr_2_add_target exe_name source_file)
  cgal_debug_message(STATUS "#   cgal_arr_2_add_target(${ARGN})")
  if(suffix)
    set(name ${exe_name}_${suffix})
  endif()
  add_executable(${name} ${source_file})
  target_link_libraries(${name} CGAL::CGAL)
  if (TARGET CGAL::CGAL_Core)
    target_link_libraries(${name} CGAL::CGAL_Core)
  endif()
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${name} )
  separate_arguments(flags UNIX_COMMAND "${TESTSUITE_CXXFLAGS}")
  target_compile_options(${name} PRIVATE ${flags})
  cgal_debug_message(STATUS "#      -> target ${name} with TESTSUITE_CXXFLAGS: ${flags}")

  if(CGAL_ENABLE_TESTING)
    cgal_add_compilation_test(${name})
  endif(CGAL_ENABLE_TESTING)

  # Add a compatibility-mode with the shell script `cgal_test_base`
  if(NOT TARGET ${exe_name})
    create_single_source_cgal_program( "${source_file}" NO_TESTING)
    if(CGAL_ENABLE_TESTING)
      cgal_add_compilation_test(${exe_name})
    endif(CGAL_ENABLE_TESTING)
  endif()
endfunction()

function(run_test name)
  # ${ARGV0} - executable name
  cgal_debug_message(STATUS "# run_test(${ARGN})")
  cgal_arr_2_add_target(${name} ${ARGN})
endfunction()

function(run_test_with_flags)
  # ${ARGV0} - executable name
  # ${ARGV1} - test substring name
  if(NOT CGAL_ENABLE_TESTING)
    return()
  endif()
  cgal_debug_message(STATUS "# run_test_with_flags(${ARGN})")
  cgal_add_test(${ARGV0}_${suffix} ${ARGV0} ${ARGV0}.${ARGV1})
endfunction()

function(run_test_alt name datafile)
  if(NOT CGAL_ENABLE_TESTING)
    return()
  endif()
  if(suffix)
    set(name ${name}_${suffix})
  endif()
  cgal_debug_message(STATUS "#     run_test_alt(${ARGN})")
  cgal_debug_message(STATUS "#       -> ./${name} ${datafile} ${ARGN}")
  string(MAKE_C_IDENTIFIER "${name}  ${ARGV4}  ${ARGV5}" test_name)
  cgal_add_test(${name}
    TEST_NAME ${test_name}
    ARGUMENTS ${datafile} ${ARGN})
endfunction()

function(run_trapped_test name datafile)
  cgal_debug_message(STATUS "#   run_trapped_test(${name} ${datafile} ${ARGN})")
  run_test_alt(${name} ${datafile} ${ARGN})
endfunction()

function(compile_and_run)
  set(name ${ARGV0})
  cgal_debug_message(STATUS "# compile_and_run(${ARGN})")
#  message("   successful compilation of ${name}")
  cgal_arr_2_add_target(${name} ${name}.cpp)
  if(CGAL_ENABLE_TESTING)
    cgal_add_test(${name})
  endif()
endfunction()

function(execute_commands_old_structure data_dir traits_type_name)
  cgal_debug_message(STATUS "# execute_commands_old_structure(data_dir=${data_dir} traits=${traits_type_name} ${ARGN})")
  # at first the tests where designed in such way that all the test input was
  # in one file, the points, the xcurves, the curves and the execution block
  # this function is used to execute the old tests, one may use it when needed
  # but you should remember that separating the input into smaller files creates
  # much more modular and comfortable test suite

  # the old structure is default, so this function executes all commands
  # except the commands that are given as arguments

  set(commands_indicator_COMPARE 1)
  set(commands_indicator_VERTEX 1)
  set(commands_indicator_IS_VERTICAL 1)
  set(commands_indicator_COMPARE_Y_AT_X 1)
  set(commands_indicator_COMPARE_Y_AT_X_LEFT 1)
  set(commands_indicator_COMPARE_Y_AT_X_RIGHT 1)
  set(commands_indicator_MAKE_X_MONOTONE 1)
  set(commands_indicator_INTERSECT 1)
  set(commands_indicator_SPLIT 1)
  set(commands_indicator_ARE_MERGEABLE 1)
  set(commands_indicator_MERGE 1)
  set(commands_indicator_ASSERTIONS 1)
  set(commands_indicator_CONSTRUCTOR 1)
  foreach(arg ${ARGN})
    set(commands_indicator_${arg} 0)
  endforeach()
  if(commands_indicator_COMPARE)
    run_trapped_test(test_traits
      data/compare.pt data/empty.zero
      data/empty.zero data/compare ${traits_type_name})
  endif()
  if(commands_indicator_VERTEX)
    run_trapped_test(test_traits
      data/${data_dir}/vertex.pt data/${data_dir}/vertex.xcv
      data/empty.zero data/${data_dir}/vertex ${traits_type_name})
  endif()
  if(commands_indicator_IS_VERTICAL)
    run_trapped_test(test_traits
      data/empty.zero data/${data_dir}/is_vertical.xcv data/empty.zero
      data/${data_dir}/is_vertical ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X)
    run_trapped_test(test_traits
      data/${data_dir}/compare_y_at_x.pt data/${data_dir}/compare_y_at_x.xcv
      data/empty.zero data/${data_dir}/compare_y_at_x ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X_LEFT)
    run_trapped_test(test_traits
      data/${data_dir}/compare_y_at_x_left.pt data/${data_dir}/compare_y_at_x_left.xcv
      data/empty.zero data/${data_dir}/compare_y_at_x_left ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X_RIGHT)
    run_trapped_test(test_traits
      data/${data_dir}/compare_y_at_x_right.pt data/${data_dir}/compare_y_at_x_right.xcv
      data/empty.zero data/${data_dir}/compare_y_at_x_right ${traits_type_name})
  endif()
  if(commands_indicator_MAKE_X_MONOTONE)
    run_trapped_test(test_traits
      data/empty.zero data/${data_dir}/make_x_monotone.xcv
      data/${data_dir}/make_x_monotone.cv data/${data_dir}/make_x_monotone ${traits_type_name})
  endif()
  if(commands_indicator_INTERSECT)
    run_trapped_test(test_traits
      data/${data_dir}/intersect.pt data/${data_dir}/intersect.xcv
      data/empty.zero data/${data_dir}/intersect ${traits_type_name})
  endif()
  if(commands_indicator_SPLIT)
    run_trapped_test(test_traits
      data/${data_dir}/split.pt data/${data_dir}/split.xcv
      data/empty.zero data/${data_dir}/split ${traits_type_name})
  endif()
  if(commands_indicator_ARE_MERGEABLE)
    run_trapped_test(test_traits
      data/empty.zero data/${data_dir}/are_mergeable.xcv
      data/empty.zero data/${data_dir}/are_mergeable ${traits_type_name})
  endif()
  if(commands_indicator_MERGE)
    run_trapped_test(test_traits
      data/empty.zero data/${data_dir}/merge.xcv
      data/empty.zero data/${data_dir}/merge ${traits_type_name})
  endif()
  if(commands_indicator_ASSERTIONS)
    run_trapped_test(test_traits
      data/${data_dir}/assertions.pt data/${data_dir}/assertions.xcv
      data/empty.zero data/${data_dir}/assertions ${traits_type_name})
  endif()
  if(commands_indicator_CONSTRUCTOR)
    run_trapped_test(test_traits
      data/empty.zero data/${data_dir}/constructor.xcv
      data/${data_dir}/constructor.cv data/${data_dir}/constructor ${traits_type_name})
  endif()
endfunction()

function(execute_commands_new_structure data_dir traits_type_name)
  cgal_debug_message(STATUS "# execute_commands_new_structure(data_dir=${data_dir} traits=${traits_type_name} ${ARGN} )")
# the new design for the tests includes separation of the test input into 4
# parts: points file, xcurves file, curves file and execution block file.
# one may reuse the input files for the various tests

# the new structure is not default, so this function executes only
# commands that are given as arguments

  set(commands_indicator_COMPARE 0)
  set(commands_indicator_VERTEX 0)
  set(commands_indicator_IS_VERTICAL 0)
  set(commands_indicator_COMPARE_X_ON_BOUNDARY 0)
  set(commands_indicator_COMPARE_X_NEAR_BOUNDARY 0)
  set(commands_indicator_COMPARE_Y_NEAR_BOUNDARY 0)
  set(commands_indicator_PARAMETER_SPACE_X 0)
  set(commands_indicator_PARAMETER_SPACE_Y 0)
  set(commands_indicator_COMPARE_Y_AT_X 0)
  set(commands_indicator_COMPARE_Y_AT_X_LEFT 0)
  set(commands_indicator_COMPARE_Y_AT_X_RIGHT 0)
  set(commands_indicator_MAKE_X_MONOTONE 0)
  set(commands_indicator_INTERSECT 0)
  set(commands_indicator_SPLIT 0)
  set(commands_indicator_ARE_MERGEABLE 0)
  set(commands_indicator_MERGE 0)
  set(commands_indicator_ASSERTIONS 0)
  set(commands_indicator_CONSTRUCTOR 0)
  set(commands_indicator_EQUAL 0)
  set(commands_indicator_PUSH_BACK 0)
  set(commands_indicator_PUSH_FRONT 0)
  set(commands_indicator_NUMBER_OF_POINTS 0)
  set(commands_indicator_COMPARE_ENDPOINTS_XY 0)
  set(commands_indicator_CONSTRUCT_OPPOSITE 0)
  set(commands_indicator_TRIM 0)
  foreach(arg ${ARGN})
    set(commands_indicator_${arg} 1)
  endforeach()
  if(commands_indicator_COMPARE)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare ${traits_type_name})
  endif()
  if(commands_indicator_VERTEX)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/vertex ${traits_type_name})
  endif()
  if(commands_indicator_IS_VERTICAL)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/is_vertical ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_X_ON_BOUNDARY)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_x_on_boundary ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_X_NEAR_BOUNDARY)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_x_near_boundary ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_NEAR_BOUNDARY)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_y_near_boundary ${traits_type_name})
  endif()
  if(commands_indicator_PARAMETER_SPACE_X)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/parameter_space_x ${traits_type_name})
  endif()
  if(commands_indicator_PARAMETER_SPACE_Y)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/parameter_space_y ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_y_at_x ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X_LEFT)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_y_at_x_left ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X_RIGHT)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_y_at_x_right ${traits_type_name})
  endif()
  if(commands_indicator_MAKE_X_MONOTONE)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/make_x_monotone ${traits_type_name})
  endif()
  if(commands_indicator_INTERSECT)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/intersect ${traits_type_name})
  endif()
  if(commands_indicator_SPLIT)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/split ${traits_type_name})
  endif()
  if(commands_indicator_ARE_MERGEABLE)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/are_mergeable ${traits_type_name})
  endif()
  if(commands_indicator_MERGE)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/merge ${traits_type_name})
  endif()
  if(commands_indicator_ASSERTIONS)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/assertions ${traits_type_name})
  endif()
  if(commands_indicator_CONSTRUCTOR)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/constructor ${traits_type_name})
  endif()
  if(commands_indicator_EQUAL)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/equal ${traits_type_name})
  endif()
  if(commands_indicator_PUSH_BACK)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/push_back ${traits_type_name})
  endif()
  if(commands_indicator_PUSH_FRONT)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/push_front ${traits_type_name})
  endif()
  if(commands_indicator_NUMBER_OF_POINTS)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/number_of_points ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_ENDPOINTS_XY)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/compare_endpoints_xy ${traits_type_name})
  endif()
  if(commands_indicator_CONSTRUCT_OPPOSITE)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/construct_opposite ${traits_type_name})
  endif()
  if(commands_indicator_TRIM)
    run_trapped_test(test_traits data/${data_dir}/points
      data/${data_dir}/xcurves data/${data_dir}/curves data/${data_dir}/trim ${traits_type_name})
  endif()
endfunction()

function(execute_commands_traits_adaptor data_dir traits_type_name)
  cgal_debug_message(STATUS "# execute_commands_traits_adaptor(data_dir=${data_dir} traits=${traits_type_name} ${ARGN} )")
# the new structure is not default, so this function executes only
# commands that are given as arguments

  set(commands_indicator_PARAMETER_SPACE_X 0)
  set(commands_indicator_PARAMETER_SPACE_Y 0)
  set(commands_indicator_COMPARE_XY 0)
  set(commands_indicator_COMPARE_X_ON_BOUNDARY 0)
  set(commands_indicator_COMPARE_X_NEAR_BOUNDARY 0)
  set(commands_indicator_COMPARE_Y_NEAR_BOUNDARY 0)
  set(commands_indicator_COMPARE_Y_AT_X_LEFT 0)
  set(commands_indicator_ARE_MERGEABLE 0)
  set(commands_indicator_MERGE 0)
  set(commands_indicator_X_ON_IDENTIFICATION 0)
  set(commands_indicator_Y_ON_IDENTIFICATION 0)
  set(commands_indicator_IS_BOUNDED 0)
  set(commands_indicator_IS_IN_X_RANGE 0)
  set(commands_indicator_COMPARE_Y_POSITION 0)
  set(commands_indicator_IS_BETWEEN_CW 0)
  set(commands_indicator_COMPARE_CW_AROUND_POINT 0)
  foreach(arg ${ARGN})
    set(commands_indicator_${arg} 1)
  endforeach()

  if(commands_indicator_PARAMETER_SPACE_X)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/parameter_space_x ${traits_type_name})
  endif()
  if(commands_indicator_PARAMETER_SPACE_Y)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/parameter_space_y ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_XY)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_xy ${traits_type_name})
  endif()

  if(commands_indicator_COMPARE_X_ON_BOUNDARY)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_x_on_boundary ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_X_NEAR_BOUNDARY)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_x_near_boundary ${traits_type_name})
  endif()

  if(commands_indicator_COMPARE_Y_NEAR_BOUNDARY)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_y_near_boundary ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_AT_X_LEFT)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_y_at_x_left ${traits_type_name})
  endif()
  if(commands_indicator_ARE_MERGEABLE)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/are_mergeable ${traits_type_name})
  endif()
  if(commands_indicator_MERGE)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/merge ${traits_type_name})
  endif()
  if(commands_indicator_X_ON_IDENTIFICATION)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/x_on_idintification ${traits_type_name})
  endif()
  if(commands_indicator_Y_ON_IDENTIFICATION)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/x_on_idintification ${traits_type_name})
  endif()
  if(commands_indicator_IS_BOUNDED)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/is_bounded ${traits_type_name})
  endif()
  if(commands_indicator_IS_IN_X_RANGE)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/is_in_x_range ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_Y_POSITION)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_y_position ${traits_type_name})
  endif()
  if(commands_indicator_IS_BETWEEN_CW)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/is_between_cw ${traits_type_name})
  endif()
  if(commands_indicator_COMPARE_CW_AROUND_POINT)
    run_trapped_test(test_traits_adaptor data/test_adaptor/${data_dir}/points
      data/test_adaptor/${data_dir}/xcurves data/test_adaptor/${data_dir}/curves
      data/test_adaptor/${data_dir}/compare_cw_around_point ${traits_type_name})
  endif()
endfunction()

#---------------------------------------------------------------------#
# traits adaptor (segments traits)
#---------------------------------------------------------------------#
function(test_segment_traits_adaptor)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits_adaptor segments "${flags}")
#  if [ -n "${SUCCESS}" ] ; then
  execute_commands_traits_adaptor( segments segments_traits_adaptor
    COMPARE_XY COMPARE_Y_POSITION COMPARE_CW_AROUND_POINT COMPARE_Y_AT_X_LEFT
    ARE_MERGEABLE MERGE IS_IN_X_RANGE IS_BETWEEN_CW)
endfunction()

#---------------------------------------------------------------------#
# traits adaptor (linear traits)
#---------------------------------------------------------------------#
function(test_linear_traits_adaptor)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags( test_traits_adaptor linear "${flags}")

  execute_commands_traits_adaptor( linear linear_traits_adaptor
    COMPARE_XY COMPARE_Y_AT_X_LEFT ARE_MERGEABLE MERGE IS_IN_X_RANGE
    COMPARE_Y_POSITION IS_BETWEEN_CW COMPARE_CW_AROUND_POINT)
endfunction()

#---------------------------------------------------------------------#
# traits adaptor (spherical arcs traits)
#---------------------------------------------------------------------#
function(test_spherical_arcs_traits_adaptor)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")

  compile_test_with_flags( test_traits_adaptor geodesic_arcs_on_sphere "${flags}")

  execute_commands_traits_adaptor( spherical_arcs spherical_arcs_traits_adaptor
    COMPARE_XY COMPARE_Y_AT_X_LEFT ARE_MERGEABLE MERGE IS_IN_X_RANGE
    COMPARE_Y_POSITION IS_BETWEEN_CW COMPARE_CW_AROUND_POINT)
endfunction()

#---------------------------------------------------------------------#
# compile and run test with traits
#---------------------------------------------------------------------#
function(compile_and_run_with_flags name type flags)
  compile_test_with_flags( ${name} ${type} "${flags}" ${ARGN})
#  if [ -n "${SUCCESS}" ] ; then
#    if [ -n "${DO_RUN}" ] ; then
  run_test_with_flags( ${name} ${type})
#    fi
#  else
#    echo "   ERROR:    not executed construction of segments" >> ${ERRORFILE}
#  fi
#  clean_tests
endfunction()

#---------------------------------------------------------------------#
# construction with segments
#---------------------------------------------------------------------#
function(test_construction_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags( test_construction segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# construction with linear curves
#---------------------------------------------------------------------#
function(test_construction_linear_curves)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(topol_traits ${PLANAR_UNBOUNDED_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")
  compile_and_run_with_flags( test_construction linear "${flags}")
endfunction()

#---------------------------------------------------------------------#
# construction with geodesic arcs on the sphere
#---------------------------------------------------------------------#
function(test_construction_spherical_arcs)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")
  compile_and_run_with_flags( test_construction geodesic_arcs_on_sphere "${flags}")
endfunction()

#---------------------------------------------------------------------#
# construction with polylines
#---------------------------------------------------------------------#
function(test_construction_polylines)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${POLYLINE_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags( test_construction polylines "${flags}")
endfunction()


#---------------------------------------------------------------------#
# overlay with segments
#---------------------------------------------------------------------#
function(test_overlay_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_overlay segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# overlay with geodesic arcs on the sphere
#---------------------------------------------------------------------#
function(test_overlay_spherical_arcs)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")
  compile_and_run_with_flags(test_overlay geodesic_arcs_on_sphere "${flags}")
endfunction()

#---------------------------------------------------------------------#
# point location with segments
#---------------------------------------------------------------------#
function(test_point_location_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_point_location segments "${flags}" segments)
endfunction()

# For backward compatibility
function(test_point_location_segments_conversion)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DCGAL_ARR_POINT_LOCATION_CONVERSION")
  compile_and_run_with_flags(test_point_location segments "${flags}" segment_conversion)
endfunction()

#---------------------------------------------------------------------#
# point location dynamic with segments
#---------------------------------------------------------------------#
function(test_point_location_dynamic_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_point_location_dynamic segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# point location with circle segments
#---------------------------------------------------------------------#
function(test_point_location_circle_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${CIRCLE_SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_point_location circle_segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# point location with linear objects
#---------------------------------------------------------------------#
function(test_point_location_linear)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_point_location linear "${flags}")
endfunction()

#---------------------------------------------------------------------#
# batchecd point location with segments
#---------------------------------------------------------------------#
function(test_batched_point_location_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_batched_point_location segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# batchecd point location with linear objects
#---------------------------------------------------------------------#
function(test_batched_point_location_linear)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_batched_point_location linear "${flags}")
endfunction()

#---------------------------------------------------------------------#
# batchecd point location with geodesic arcs on the sphere
#---------------------------------------------------------------------#
function(test_batched_point_location_spherical_arcs)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")
  compile_and_run_with_flags(test_batched_point_location geodesic_arcs_on_sphere "${flags}")
endfunction()

#---------------------------------------------------------------------#
# vertical decomposition with segments
#---------------------------------------------------------------------#
function(test_vertical_decomposition_segments)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_vertical_decomposition segments "${flags}")
endfunction()

#---------------------------------------------------------------------#
# vertical decomposition with linear objects
#---------------------------------------------------------------------#
function(test_vertical_decomposition_linear)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")
  compile_and_run_with_flags(test_vertical_decomposition linear "${flags}")
endfunction()

#---------------------------------------------------------------------#
# vertical decomposition with geodesic arcs on the sphere
#---------------------------------------------------------------------#
function(test_vertical_decomposition_spherical_arcs)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")
  compile_and_run_with_flags(test_vertical_decomposition geodesic_arcs_on_sphere "${flags}")
endfunction()

#---------------------------------------------------------------------#
# segment traits
#---------------------------------------------------------------------#
function(test_segment_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits segments "${flags}")

  execute_commands_old_structure( segments segment_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT CONSTRUCTOR
    COMPARE_Y_AT_X_RIGHT ARE_MERGEABLE)

  execute_commands_new_structure( segments segment_traits
    IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT ARE_MERGEABLE)

  run_trapped_test( test_traits
    data/segments/vertex.pt data/segments/xcurves
    data/empty.zero data/segments/vertex segment_traits)
endfunction()

#---------------------------------------------------------------------#
# non-caching segment traits
#---------------------------------------------------------------------#
function(test_non_caching_segment_traits)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${NON_CACHING_SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits non_caching_segments "${flags}")

  execute_commands_old_structure(segments non_caching_segment_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT CONSTRUCTOR
    COMPARE_Y_AT_X_RIGHT ARE_MERGEABLE ASSERTIONS)

  execute_commands_new_structure(segments segment_traits
    IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT)

  run_trapped_test(test_traits
    data/segments/vertex.pt data/segments/xcurves
    data/empty.zero data/segments/vertex non_caching_segment_traits)
endfunction()

#---------------------------------------------------------------------#
# polycurve conic traits
#---------------------------------------------------------------------#
function(test_polycurve_conic_traits)
#  echo polycurve test starting
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_polycurve_conic_traits requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_EXPR_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${POLYCURVE_CONIC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits conic_polycurve "${flags}")

  # The input arguments for the execute_commands_new_structure,
  # 1. polycurve_conics is the directory name in "data"
  # 2. polycurve_conic_traits is a string
  # Execute_command_new_structure will only run the test on functors provided as the third, fourth and so on arguments.
  # To see how the input data directory should be structured for each functor, check the execute_commands_new_structure function in this file.
  execute_commands_new_structure(polycurves_conics polycurve_conic_traits
    COMPARE_Y_AT_X
    INTERSECT
    EQUAL
    IS_VERTICAL
    SPLIT
    ARE_MERGEABLE
    COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT
    MAKE_X_MONOTONE
    PUSH_BACK
    PUSH_FRONT
    NUMBER_OF_POINTS
    VERTEX
    CONSTRUCT_OPPOSITE
    MERGE
    COMPARE_ENDPOINTS_XY
    TRIM)

endfunction()

#---------------------------------------------------------------------#
# polycurve arc traits
#---------------------------------------------------------------------#
function(test_polycurve_circular_arc_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits circular_arc_polycurve "${flags}")

  execute_commands_new_structure(polycurves_circular_arcs polycurve_circular_arc_traits
    COMPARE_Y_AT_X
    EQUAL
    IS_VERTICAL
    SPLIT
    ARE_MERGEABLE
    COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT
    MAKE_X_MONOTONE
    PUSH_BACK
    PUSH_FRONT
    NUMBER_OF_POINTS
    VERTEX
    CONSTRUCT_OPPOSITE
    MERGE
    COMPARE_ENDPOINTS_XY
    INTERSECT)
endfunction()

#---------------------------------------------------------------------#
# polycurve bezier traits
#---------------------------------------------------------------------#
function(test_polycurve_bezier_traits)
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_polycurve_bezier_traits requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_EXPR_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${POLYCURVE_BEZIER_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits bezier_polycurve "${flags}")

  execute_commands_new_structure(polycurves_bezier test_polycurve_bezier_traits
    MERGE
    EQUAL
    IS_VERTICAL
    NUMBER_OF_POINTS
    PUSH_BACK
    PUSH_FRONT
    VERTEX
    ARE_MERGEABLE
    COMPARE_ENDPOINTS_XY
    # TODO (add data for these tests)
    # COMPARE_Y_AT_X
    # SPLIT
    # COMPARE_Y_AT_X_LEFT
    # COMPARE_Y_AT_X_RIGHT
    # MAKE_X_MONOTONE
    # CONSTRUCT_OPPOSITE

    # INTERSECT
    )
endfunction()

#---------------------------------------------------------------------#
# polyline traits
#---------------------------------------------------------------------#
function(test_polyline_traits)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${POLYLINE_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits test_polylines "${flags}")

  execute_commands_old_structure(polylines polyline_traits
    CONSTRUCTOR COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT ARE_MERGEABLE)
endfunction()

#---------------------------------------------------------------------#
# non-caching polyline traits
#---------------------------------------------------------------------#
function(test_non_caching_polyline_traits)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${NON_CACHING_POLYLINE_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits non_caching_polylines "${flags}")

  execute_commands_old_structure(polylines non_caching_polyline_traits
    CONSTRUCTOR COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT ARE_MERGEABLE)
endfunction()

#---------------------------------------------------------------------#
# linear traits
#---------------------------------------------------------------------#
function(test_linear_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINEAR_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits linear "${flags}")

  execute_commands_old_structure(linear/segments linear_traits.segments
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT CONSTRUCTOR ARE_MERGEABLE)

  execute_commands_new_structure(linear/segments linear_traits.segments
    IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT)

  run_trapped_test(test_traits
    data/linear/segments/vertex.pt data/linear/segments/xcurves
    data/empty.zero data/linear/segments/vertex linear_traits.segments)

  execute_commands_old_structure(linear/rays linear_traits.rays
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT CONSTRUCTOR ARE_MERGEABLE)

  execute_commands_new_structure(linear/rays linear_traits.rays
    IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT)

  run_trapped_test(test_traits
    data/linear/rays/vertex.pt data/linear/rays/xcurves
    data/empty.zero data/linear/rays/vertex linear_traits.rays)

  execute_commands_new_structure(linear/lines linear_traits.lines
    IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT INTERSECT
    SPLIT MERGE
    PARAMETER_SPACE_X PARAMETER_SPACE_Y
    COMPARE_X_ON_BOUNDARY COMPARE_X_NEAR_BOUNDARY COMPARE_Y_NEAR_BOUNDARY)
endfunction()

#---------------------------------------------------------------------#
# conic traits
#---------------------------------------------------------------------#
function(test_conic_traits)
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_conic_traits requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_EXPR_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${CORE_CONIC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits conics "${flags}")

  execute_commands_old_structure(conics conic_traits
    INTERSECT SPLIT MERGE COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT ARE_MERGEABLE)

  execute_commands_new_structure(conics conic_traits
    INTERSECT SPLIT MERGE)

  run_trapped_test(test_traits
    data/conics/compare.pt data/empty.zero
    data/empty.zero data/conics/compare conic_traits)
endfunction()

#---------------------------------------------------------------------#
# "line arcs" (segments) only
#---------------------------------------------------------------------#
function(test_line_arc_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${LINE_ARC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits line_arcs "${flags}")

  execute_commands_old_structure(circular_lines line_arc_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    ASSERTIONS COMPARE_Y_AT_X_RIGHT MERGE ARE_MERGEABLE)

  execute_commands_new_structure(circular_lines line_arc_traits
    IS_VERTICAL COMPARE_Y_AT_X)

  run_trapped_test(test_traits
    data/circular_lines/compare.pt data/empty.zero
    data/empty.zero data/circular_lines/compare line_arc_traits)

  run_trapped_test(test_traits
    data/circular_lines/vertex.pt data/circular_lines/xcurves
    data/empty.zero data/circular_lines/vertex line_arc_traits)
endfunction()

#---------------------------------------------------------------------#
# circular arcs only
#---------------------------------------------------------------------#
function(test_circular_arc_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${CIRCULAR_ARC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits circular_arcs "${flags}")

  execute_commands_old_structure(circular_arcs circular_arc_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    ASSERTIONS COMPARE_Y_AT_X_RIGHT MERGE ARE_MERGEABLE)

  execute_commands_new_structure(circular_arcs circular_arc_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X)
endfunction()

#---------------------------------------------------------------------#
# circular and line arcs
#---------------------------------------------------------------------#
function(test_circular_line_arc_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${CIRCULAR_LINE_ARC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits circular_line_arcs "${flags}")

  execute_commands_old_structure(circular_line_arcs circular_line_arc_traits
    VERTEX IS_VERTICAL CONSTRUCTOR COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    ASSERTIONS COMPARE_Y_AT_X_RIGHT MERGE ARE_MERGEABLE)

  execute_commands_new_structure(circular_line_arcs circular_line_arc_traits
    IS_VERTICAL COMPARE_Y_AT_X)

  run_trapped_test(test_traits
    data/circular_line_arcs/vertex.pt data/circular_line_arcs/xcurves
    data/empty.zero data/circular_line_arcs/vertex circular_line_arc_traits)
endfunction()

#---------------------------------------------------------------------#
# circle segment traits
#---------------------------------------------------------------------#
function(test_circle_segments_traits)
  set(nt ${QUOTIENT_MP_FLOAT_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${CIRCLE_SEGMENT_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits circle_segments "${flags}")

  execute_commands_old_structure(circle_segments circle_segments_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT
    COMPARE_Y_AT_X_RIGHT CONSTRUCTOR ARE_MERGEABLE)

  run_trapped_test(test_traits
    data/circle_segments/points data/circle_segments/xcurves.8
    data/empty.zero data/circle_segments/vertex circle_segments_traits)
  run_trapped_test(test_traits
    data/empty.zero data/circle_segments/xcurves.8
    data/empty.zero data/circle_segments/is_vertical circle_segments_traits)
  run_trapped_test(test_traits
    data/circle_segments/points data/circle_segments/xcurves.8
    data/empty.zero data/circle_segments/compare_y_at_x circle_segments_traits)
  run_trapped_test(test_traits
    data/circle_segments/points data/circle_segments/xcurves.16
    data/empty.zero data/circle_segments/compare_y_at_x_left circle_segments_traits)
  run_trapped_test(test_traits
    data/circle_segments/points data/circle_segments/xcurves.16
    data/empty.zero data/circle_segments/compare_y_at_x_right circle_segments_traits)
  run_trapped_test(test_traits
    data/empty.zero data/circle_segments/constructor.xcv
    data/empty.zero data/circle_segments/constructor circle_segments_traits)
endfunction()

#---------------------------------------------------------------------#
# bezier traits
#---------------------------------------------------------------------#
function(test_bezier_traits)
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_bezier_traits requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_EXPR_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${BEZIER_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits Bezier "${flags}")

  execute_commands_old_structure(bezier bezier_traits
    COMPARE_Y_AT_X_LEFT COMPARE_Y_AT_X_RIGHT SPLIT
    CONSTRUCTOR ASSERTIONS ARE_MERGEABLE)
endfunction()

#---------------------------------------------------------------------#
# spherical arc traits
#---------------------------------------------------------------------#
function(test_spherical_arc_traits)
  set(nt ${CGAL_GMPQ_NT})
  set(kernel ${CARTESIAN_KERNEL})
  set(geom_traits ${GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS})
  set(topol_traits ${SPHERICAL_TOPOL_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits} -DTEST_TOPOL_TRAITS=${topol_traits}")

  compile_test_with_flags(test_traits geodesic_arcs_on_sphere "${flags}")

  execute_commands_old_structure(spherical_arcs spherical_arc_traits
    COMPARE_Y_AT_X_LEFT COMPARE_Y_AT_X_RIGHT INTERSECT
    CONSTRUCTOR
    COMPARE MAKE_X_MONOTONE SPLIT MERGE ASSERTIONS ARE_MERGEABLE)

  execute_commands_new_structure(spherical_arcs spherical_arc_traits
    INTERSECT
    COMPARE_X_ON_BOUNDARY COMPARE_X_NEAR_BOUNDARY
    COMPARE_Y_NEAR_BOUNDARY)

  run_trapped_test(test_traits
    data/spherical_arcs/compare.pt data/spherical_arcs/compare.xcv
    data/empty.zero data/spherical_arcs/compare spherical_arc_traits)
endfunction()

#---------------------------------------------------------------------#
# rational arc traits
#---------------------------------------------------------------------#
function(test_rational_arc_traits)
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_rational_arc_traits requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_INT_NT})
  set(kernel ${UNIVARIATE_ALGEBRAIC_KERNEL})
  set(geom_traits ${RATIONAL_ARC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits rational_arcs "${flags}")

  run_trapped_test(test_traits
    data/compare.pt data/empty.zero
    data/empty.zero data/compare rational_arc_traits)

  execute_commands_new_structure(rational_arcs rational_arc_traits
    VERTEX IS_VERTICAL COMPARE_Y_AT_X COMPARE_Y_AT_X_LEFT SPLIT MERGE
    COMPARE_X_ON_BOUNDARY COMPARE_X_NEAR_BOUNDARY COMPARE_Y_NEAR_BOUNDARY)
endfunction()

#---------------------------------------------------------------------#
# algebraic traits with GMP/MPFI
#---------------------------------------------------------------------#
function(test_algebraic_traits_gmp)
  #TODO: Adapt
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_traits_algebraic_traits_gmp requires GMP and will not be executed")
    return()
  endif()
  set(nt ${CGAL_GMPZ_NT})
  set(kernel ${UNIVARIATE_ALGEBRAIC_KERNEL})
  set(geom_traits ${ALGEBRAIC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits algebraic "${flags}" algebraic_traits_gmp)

  execute_commands_new_structure(algebraic algebraic_traits_gmp
    COMPARE COMPARE_Y_AT_X COMPARE_Y_AT_X_RIGHT COMPARE_Y_AT_X_LEFT
    MAKE_X_MONOTONE IS_VERTICAL VERTEX SPLIT MERGE INTERSECT
    PARAMETER_SPACE_X PARAMETER_SPACE_Y)
endfunction()

#---------------------------------------------------------------------#
# algebraic traits with LEDA
#---------------------------------------------------------------------#
function(test_algebraic_traits_leda)
  #TODO: Adapt

  set(nt ${LEDA_INT_NT})
  set(kernel ${UNIVARIATE_ALGEBRAIC_KERNEL})
  set(geom_traits ${ALGEBRAIC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits algebraic "${flags}" algebraic_traits_leda)

  execute_commands_new_structure(algebraic algebraic_traits_leda
    COMPARE COMPARE_Y_AT_X COMPARE_Y_AT_X_RIGHT COMPARE_Y_AT_X_LEFT
    MAKE_X_MONOTONE IS_VERTICAL VERTEX SPLIT MERGE INTERSECT
    PARAMETER_SPACE_X PARAMETER_SPACE_Y)
endfunction()


#---------------------------------------------------------------------#
# algebraic traits with CORE
#---------------------------------------------------------------------#
function(test_algebraic_traits_core)
  #TODO: Adapt
  if(CGAL_DISABLE_GMP)
    MESSAGE(STATUS "test_algebraic_traits_core requires CORE and will not be executed")
    return()
  endif()
  set(nt ${CORE_INT_NT})
  set(kernel ${UNIVARIATE_ALGEBRAIC_KERNEL})
  set(geom_traits ${ALGEBRAIC_GEOM_TRAITS})
  set(flags "-DTEST_NT=${nt} -DTEST_KERNEL=${kernel} -DTEST_GEOM_TRAITS=${geom_traits}")

  compile_test_with_flags(test_traits algebraic "${flags}" algebraic_traits_core)

  execute_commands_new_structure(algebraic algebraic_traits_core
    COMPARE COMPARE_Y_AT_X COMPARE_Y_AT_X_RIGHT COMPARE_Y_AT_X_LEFT
    MAKE_X_MONOTONE IS_VERTICAL VERTEX SPLIT MERGE INTERSECT
    PARAMETER_SPACE_X PARAMETER_SPACE_Y)
endfunction()


configure("")
compile_and_run(construction_test_suite_generator)

test_segment_traits()
test_non_caching_segment_traits()
test_polyline_traits()
test_polycurve_conic_traits()
test_polycurve_circular_arc_traits()
test_polycurve_bezier_traits()
test_non_caching_polyline_traits()
test_linear_traits()
test_conic_traits()

test_line_arc_traits() 		# "line arcs" (segments) only
test_circular_arc_traits() 	# circular arcs only
test_circular_line_arc_traits() 	# for both

test_circle_segments_traits()
test_bezier_traits()

test_spherical_arc_traits()

test_rational_arc_traits()

test_algebraic_traits_core()
test_algebraic_traits_gmp()
test_algebraic_traits_leda()

compile_and_run(test_data_traits)

compile_and_run(test_insertion)
compile_and_run(test_unbounded_rational_insertion)
compile_and_run(test_unbounded_rational_direct_insertion)
compile_and_run(test_rational_function_traits_2)
compile_and_run(test_iso_verts)

compile_and_run(test_vert_ray_shoot_vert_segments)

test_construction_segments()
test_construction_linear_curves()
test_construction_spherical_arcs()
test_construction_polylines()

test_overlay_segments()
test_overlay_spherical_arcs()

test_point_location_segments()
test_point_location_segments_conversion()
test_point_location_circle_segments()
test_point_location_linear()

test_point_location_dynamic_segments()

test_batched_point_location_segments()
test_batched_point_location_linear()
test_batched_point_location_spherical_arcs()

test_vertical_decomposition_segments()
test_vertical_decomposition_linear()
# test_vertical_decomposition_spherical_arcs

compile_and_run(test_dual)
compile_and_run(test_do_intersect)
compile_and_run(test_zone)

compile_and_run(test_observer)
compile_and_run(test_do_equal)

test_segment_traits_adaptor()
test_linear_traits_adaptor()
test_spherical_arcs_traits_adaptor()

compile_and_run(test_removal)
compile_and_run(test_unbounded_removal)
compile_and_run(test_spherical_removal)

compile_and_run(test_io)

compile_and_run(test_sgm)

compile_and_run(test_polycurve_intersection)
if(CGAL_DISABLE_GMP)
  get_directory_property(LIST_OF_TESTS TESTS)
  foreach(_test ${LIST_OF_TESTS})
    set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT CGAL_DISABLE_GMP=1)
  endforeach()
endif()
