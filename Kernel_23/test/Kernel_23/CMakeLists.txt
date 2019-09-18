# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Kernel_23_Tests )

set(cppfiles
  Cartesian.cpp
  Dimension.cpp
  Exact_predicates_exact_constructions_kernel.cpp
  Filtered_cartesian.cpp
  Filtered_homogeneous.cpp
  Homogeneous.cpp
  issue_129.cpp
  issue_3301.cpp
  Kernel_checker.cpp
  Lazy_kernel.cpp
  overload_bug.cpp
  Simple_cartesian.cpp
  Simple_homogeneous.cpp
  test_all_linear_intersections.cpp
  test_bbox.cpp
  test_Has_conversion.cpp
  test_hash_functions.cpp
  test_kernel__.cpp
  test_projection_traits.cpp
  test_Projection_traits_xy_3_Intersect_2.cpp
  )

find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include_directories (BEFORE "include")

  # create a target per cppfile
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()
  find_package(benchmark)
  if(TARGET benchmark::benchmark)
    foreach(bench
        bench_Epeck
        bench_Epeck-no-atomic
        bench_Epeck-with-atomic-intervals
        bench_Epeck-with-mutex
        bench_Epeck-with-shared-timed-mutex
        bench_Epeck-with-shared-mutex
        )
      add_executable(${bench} bench_Epeck.cpp)
      target_link_libraries(${bench} PRIVATE CGAL::CGAL benchmark::benchmark)
      cgal_add_test(${bench})
    endforeach()
    target_compile_definitions(bench_Epeck-with-mutex              PRIVATE CGAL_LAZY_USE_MUTEX)
    target_compile_definitions(bench_Epeck-with-shared-mutex       PRIVATE CGAL_LAZY_USE_SHARED_MUTEX)
    target_compile_definitions(bench_Epeck-with-shared-timed-mutex PRIVATE CGAL_LAZY_USE_SHARED_TIMED_MUTEX)
    target_compile_definitions(bench_Epeck-with-atomic-intervals   PRIVATE CGAL_LAZY_USE_ATOMIC_APPROX)

    target_compile_features   (bench_Epeck-with-shared-mutex       PRIVATE cxx_std_17)
  else()
    message(FATAL_ERROR "Cannot find Google benchmark")
  endif()
else()
  
    message(STATUS "This program requires the CGAL and CGAL_Core libraries, and will not be compiled.")
  
endif()

